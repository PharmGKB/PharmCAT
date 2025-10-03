package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.BaseConfig;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.Constants.isLowestFunctionGene;


/**
 * This is the main entry point for matching {@link NamedAllele}s.
 *
 * @author Mark Woon
 */
public class NamedAlleleMatcher {
  public static final String VERSION = "2.0.0";
  // CHANGES TO THIS LIST OF GENES SHOULD ALSO BE REFLECTED IN DOCUMENTATION IN NamedAlleleMatcher-101.md
  public static final List<String> TREAT_UNDOCUMENTED_VARIATIONS_AS_REFERENCE = List.of(
      "CACNA1S",
      "DPYD",
      "G6PD",
      "NUDT15",
      "RYR1",
      "TPMT"
  );

  private final Env m_env;
  private final DefinitionReader m_definitionReader;
  private final boolean m_findCombinations;
  private final boolean m_topCandidateOnly;
  private final boolean m_callCyp2d6;
  private boolean m_printWarnings;


  /**
   * Default constructor.
   * This will only call the top candidate(s).  Will not find combinations or call CYP2D6.
   */
  public NamedAlleleMatcher(Env env, DefinitionReader definitionReader) {
    this(env, definitionReader, false, false, false);
  }

  /**
   * Constructor.
   *
   * @param topCandidateOnly true if only top candidate(s) should be called, false to call all possible candidates
   * @param callCyp2d6 true if CYP2D6 should be called
   */
  public NamedAlleleMatcher(Env env, DefinitionReader definitionReader, boolean findCombinations,
      boolean topCandidateOnly, boolean callCyp2d6) {
    Preconditions.checkNotNull(env);
    Preconditions.checkNotNull(definitionReader);
    m_env = env;
    m_definitionReader = definitionReader;
    m_findCombinations = findCombinations;
    m_topCandidateOnly = topCandidateOnly;
    m_callCyp2d6 = callCyp2d6;
  }


  public NamedAlleleMatcher printWarnings() {
    m_printWarnings = true;
    return this;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("vcf", "matcher-vcf", "input VCF file", true, "vcf")
          .addOption("a", "all-results", "return all possible diplotypes, not just top hits")
          // optional data
          .addOption("d", "definitions-dir", "directory containing named allele definitions (JSON files)", false, "dir")
          // outputs
          .addOption("o", "output-dir", "directory to output to (optional, default is input file directory)", false, "o")
          .addOption("html", "save-html", "save matcher results as HTML")
          // research
          .addOption("research", "research-mode", "enable research mode")
          .addOption("cyp2d6", "research-cyp2d6", "call CYP2D6 (must also use research mode)")
          .addOption("combinations", "research-combinations", "find combinations and partial alleles (must also use research mode)")
          ;

      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path vcfFile = cliHelper.getValidFile("vcf", false);
      Path definitionDir;
      if (cliHelper.hasOption("d")) {
        definitionDir = cliHelper.getValidDirectory("d", false);
      } else {
        definitionDir = DataManager.DEFAULT_DEFINITION_DIR;
      }

      DefinitionReader definitionReader = new DefinitionReader(definitionDir);
      if (definitionReader.getGenes().isEmpty()) {
        System.out.println("Did not find any allele definitions at " + definitionDir);
        System.exit(1);
      }

      boolean topCandidateOnly = !cliHelper.hasOption("a");
      boolean callCyp2d6 = false;
      boolean findCombinations = false;
      if (cliHelper.hasOption("research")) {
        callCyp2d6 = cliHelper.hasOption("cyp2d6");
        findCombinations = cliHelper.hasOption("combinations");
      }
      NamedAlleleMatcher namedAlleleMatcher =
          new NamedAlleleMatcher(new Env(), definitionReader, findCombinations, topCandidateOnly, callCyp2d6)
              .printWarnings();
      Result result = namedAlleleMatcher.call(new VcfFile(vcfFile), null, null);

      Path jsonFile = CliUtils.getOutputFile(cliHelper, vcfFile, "json", BaseConfig.MATCHER_SUFFIX + ".json");
      ResultSerializer resultSerializer = new ResultSerializer();
      resultSerializer.toJson(result, jsonFile);
      System.out.println("Saved JSON results to " + jsonFile);
      if (cliHelper.hasOption("html")) {
        Path htmlFile = CliUtils.getOutputFile(cliHelper, vcfFile, "html", BaseConfig.MATCHER_SUFFIX + ".html");
        resultSerializer.toHtml(result, htmlFile);
        System.out.println("Saved HTML results to " + htmlFile);
      }

    } catch (Exception ex) {
      //noinspection CallToPrintStackTrace
      ex.printStackTrace();
    }
  }


  public void saveResults(Result result, @Nullable Path jsonFile, @Nullable Path htmlFile) throws IOException {
    ResultSerializer resultSerializer = new ResultSerializer();
    if (jsonFile != null) {
      resultSerializer.toJson(result, jsonFile);
    }
    if (htmlFile != null) {
      resultSerializer.toHtml(result, htmlFile);
    }
  }


  /**
   * Calls diplotypes for the given VCF file for all genes for which a definition exists.
   */
  public Result call(VcfFile vcfFile, @Nullable String sampleId) throws IOException {
    return call(vcfFile, sampleId, null);
  }

  /**
   * Calls diplotypes for the given VCF file for all genes for which a definition exists.
   */
  public Result call(VcfFile vcfFile, @Nullable String sampleId, @Nullable Path sampleMetadataFile) throws IOException {
    VcfReader vcfReader = vcfFile.getReader(m_definitionReader, sampleId, m_findCombinations);
    SortedMap<String, SampleAllele> alleleMap = vcfReader.getAlleleMap();
    ResultBuilder resultBuilder = new ResultBuilder(m_definitionReader, m_topCandidateOnly, m_findCombinations, m_callCyp2d6)
        .forFile(vcfFile, vcfReader.getWarnings().asMap(), vcfReader.getSampleId(), sampleMetadataFile);

    if (m_printWarnings) {
      vcfReader.getWarnings().keySet()
          .forEach(key -> {
            System.out.println(key);
            vcfReader.getWarnings().get(key)
                .forEach(msg -> System.out.println("\t" + msg));
          });
    }
    // call haplotypes
    for (String gene : m_definitionReader.getGenes()) {
      if (!m_callCyp2d6 && gene.equals("CYP2D6")) {
        continue;
      }
      if (isLowestFunctionGene(gene)) {
        callLowestFunctionGene(vcfReader.getSampleId(), gene, alleleMap, resultBuilder);
      } else {
        callAssumingReference(vcfReader.getSampleId(), alleleMap, gene, resultBuilder);
      }
    }
    return resultBuilder.build(m_env);
  }

  private boolean hasRequiredPositions(MatchData data, ResultBuilder resultBuilder) {
    if (!data.getMissingRequiredPositions().isEmpty()) {
      StringBuilder builder = new StringBuilder("Cannot call ")
          .append(data.getGene())
          .append(" - missing required variant");
      if (data.getMissingRequiredPositions().size() > 1) {
        builder.append("s");
      }
      builder.append(" (")
          .append(String.join(", ", data.getMissingRequiredPositions()))
          .append(")");
      resultBuilder.noCall(data.getGene(), data, List.of(new MessageAnnotation(MessageAnnotation.TYPE_NOTE,
          "missing-required-position", builder.toString())));
      return false;
    }
    // AMP 1 required positions check is only a warning and is handled by ResultBuilder
    return true;
  }

  /**
   * Call standard gene haplotypes.
   * Missing alleles in {@link NamedAllele}s should be treated as reference.
   * This determines how to interpret PharmCAT named allele definitions, not sample data.
   */
  private void callAssumingReference(String sampleId, SortedMap<String, SampleAllele> alleleMap, String gene,
      ResultBuilder resultBuilder) {

    MatchData data = initializeCallData(sampleId, alleleMap, gene, true, false);
    if (data.getNumSampleAlleles() == 0) {
      resultBuilder.noCall(gene, data);
      return;
    }
    if (!hasRequiredPositions(data, resultBuilder)) {
      return;
    }

    if (data.hasPartialMissingAlleles()) {
      if (m_findCombinations) {
        callCombination(sampleId, alleleMap, gene, resultBuilder);
      } else {
        resultBuilder.noCall(gene, data);
      }
      return;
    }

    SortedSet<DiplotypeMatch> matches = new DiplotypeMatcher(m_env, data)
        .compute(false, m_topCandidateOnly);
    if (matches.isEmpty()) {
      if (!m_findCombinations) {
        resultBuilder.diplotypes(gene, data, matches);
      } else {
        callCombination(sampleId, alleleMap, gene, resultBuilder);
      }
      return;
    }
    resultBuilder.diplotypes(gene, data, matches);
  }

  private void callCombination(String sampleId, SortedMap<String, SampleAllele> alleleMap, String gene,
      ResultBuilder resultBuilder) {

    MatchData data = initializeCallData(sampleId, alleleMap, gene, false, true);
    if (data.getNumSampleAlleles() == 0) {
      resultBuilder.noCall(gene, data);
      return;
    }

    SortedSet<DiplotypeMatch> matches = new DiplotypeMatcher(m_env, data)
        .compute(true, false);
    resultBuilder.diplotypes(gene, data, matches);
  }


  /**
   * Calls diplotype based on the lowest-function algorithm:
   * <ul>
   *   <li>Only attempt diplotype match if effectively phased</li>
   *   <li>
   *     If no exact match or not effectively phased, only look for potential haplotypes
   *     (try to match all permutations to any potential haplotype without assuming reference)
   *   </li>
   * </ul>
   *
   */
  private void callLowestFunctionGene(String sampleId, String gene, SortedMap<String, SampleAllele> alleleMap,
      ResultBuilder resultBuilder) {

    MatchData origData = initializeCallData(sampleId, alleleMap, gene, true, false);
    if (origData.getNumSampleAlleles() == 0) {
      resultBuilder.noCall(gene, origData);
      return;
    }
    if (!hasRequiredPositions(origData, resultBuilder)) {
      return;
    }

    DpydHapB3Matcher dpydHapB3Matcher = null;
    MatchData workingData = origData;
    if (gene.equals("DPYD")) {
      dpydHapB3Matcher = new DpydHapB3Matcher(m_env, alleleMap, origData);
      if (dpydHapB3Matcher.hasHapB3Variants()) {
        workingData = initializeDpydCallData(sampleId, alleleMap, true, false);
      }
    }

    // try for diplotypes if effectively phased
    if (origData.isEffectivelyPhased()) {
      // look for exact matches (use topCandidateOnly = false because looking for exact match)
      SortedSet<DiplotypeMatch> diplotypeMatches = new DiplotypeMatcher(m_env, workingData)
          .compute(false, false);
      if (dpydHapB3Matcher != null && dpydHapB3Matcher.hasHapB3Variants()) {
        // must add HapB3 call
        if (!diplotypeMatches.isEmpty()) {
          // has matches, add HapB3 call
          MatchData mergerData = initializeCallData(sampleId, alleleMap, gene, false, false);
          SortedSet<DiplotypeMatch> mergedMatches = dpydHapB3Matcher.mergePhasedHapB3Call(mergerData, diplotypeMatches);
          resultBuilder.diplotypes(gene, mergerData, mergedMatches, dpydHapB3Matcher.getWarnings());
          return;
        } else if (!dpydHapB3Matcher.hasNonHapB3Variants()) {
          // no matches (so reference everywhere else) - only has HapB3
          MatchData mergerData = initializeCallData(sampleId, alleleMap, gene, false, false);
          SortedSet<DiplotypeMatch> mergedMatches = dpydHapB3Matcher.addPhasedHapB3CallToRef(mergerData);
          resultBuilder.diplotypes(gene, mergerData, mergedMatches, dpydHapB3Matcher.getWarnings());
          return;
        }
      }
      // default behavior (including DPYD when there are no HapB3 variants)
      if (!diplotypeMatches.isEmpty()) {
        resultBuilder.diplotypes(gene, workingData, diplotypeMatches);
        return;
      }
    }

    // try for diplotypes with combinations
    if (origData.isEffectivelyPhased() || origData.isUsingPhaseSets()) {
      MatchData comboData;
      if (dpydHapB3Matcher != null && dpydHapB3Matcher.hasHapB3Variants()) {
        comboData = initializeDpydCallData(sampleId, alleleMap, false, true);
      } else {
        comboData = initializeCallData(sampleId, alleleMap, gene, false, true);
      }
      // look for combinations
      SortedSet<DiplotypeMatch> diplotypeMatches = new DiplotypeMatcher(m_env, comboData)
          .compute(true, false);
      if (dpydHapB3Matcher != null && dpydHapB3Matcher.hasHapB3Variants()) {
        // must add HapB3 call
        if (!diplotypeMatches.isEmpty()) {
          MatchData mergerData = initializeCallData(sampleId, alleleMap, gene, false, false);
          diplotypeMatches = dpydHapB3Matcher.fixPartials(mergerData, diplotypeMatches);
          SortedSet<DiplotypeMatch> mergedMatches = dpydHapB3Matcher.mergePhasedHapB3Call(mergerData, diplotypeMatches);
          resultBuilder.diplotypes(gene, mergerData, mergedMatches, dpydHapB3Matcher.getWarnings());
          return;
        } else if (!dpydHapB3Matcher.hasNonHapB3Variants()) {
          // no matches (so reference everywhere else) - only has HapB3
          MatchData mergerData = initializeCallData(sampleId, alleleMap, gene, false, false);
          SortedSet<DiplotypeMatch> mergedMatches = dpydHapB3Matcher.addPhasedHapB3CallToRef(mergerData);
          resultBuilder.diplotypes(gene, mergerData, mergedMatches, dpydHapB3Matcher.getWarnings());
          return;
        }
      }
      // default behavior (including DPYD when there are no HapB3 variants)
      if (!diplotypeMatches.isEmpty()) {
        resultBuilder.diplotypes(gene, comboData, diplotypeMatches);
        return;
      }
    }

    // try combinations
    MatchData comboData;
    if (dpydHapB3Matcher != null) {
      comboData = initializeDpydCallData(sampleId, alleleMap, false, true);
    } else {
      comboData = initializeCallData(sampleId, alleleMap, gene, false, true);
    }
    SortedSet<DiplotypeMatch> comboDipMatches = new DiplotypeMatcher(m_env, comboData)
        .compute(true, false, false);
    if (!comboDipMatches.isEmpty()) {
      if (origData.isEffectivelyPhased() || comboData.isUsingPhaseSets()) {
        if (dpydHapB3Matcher != null && dpydHapB3Matcher.hasHapB3Variants()) {
          MatchData mergerData = initializeCallData(sampleId, alleleMap, gene, false, true);
          SortedSet<DiplotypeMatch> mergedMatches = dpydHapB3Matcher.mergePhasedHapB3Call(mergerData, comboDipMatches);
          resultBuilder.diplotypes(gene, mergerData, mergedMatches, dpydHapB3Matcher.getWarnings());
        } else {
          resultBuilder.diplotypes(gene, comboData, comboDipMatches);
        }
        return;
      } else {
        // any unphased diplotype matches would be handled here
        // TODO(markwoon): consider reporting all diplotypes, and letting reporter pick lowest function
        // TODO(markwoon): this is a more reliable way of dealing with possible strand combinations
        // the following commented out code would merge HapB3 calls into the diplotype results,
        // but is not compatible with how DpydHapB3Matcher.callHapBeHaplotypes() works
//      } else if (dpydHapB3Matcher != null && dpydHapB3Matcher.hasHapB3Variants()) {
//        comboData = initializeCallData(sampleId, alleleMap, gene, false, true);
//        comboDipMatches = dpydHapB3Matcher.mergePhasedHapB3Call(comboData, comboDipMatches);
      }
    }

    List<HaplotypeMatch> hapMatches = callHaplotypesForLowestFunctionGene(comboData, comboDipMatches, dpydHapB3Matcher);
    if (dpydHapB3Matcher != null) {
      resultBuilder.haplotypes(gene, origData, hapMatches, dpydHapB3Matcher.getWarnings());
    } else {
      resultBuilder.haplotypes(gene, origData, hapMatches);
    }
  }


  private List<HaplotypeMatch> callHaplotypesForLowestFunctionGene(MatchData comboData,
      SortedSet<DiplotypeMatch> matches, @Nullable DpydHapB3Matcher dpydHapB3Matcher) {

    if (dpydHapB3Matcher != null) {
      // must call regardless of whether we have HapB3Variants because it initializes variables used later
      dpydHapB3Matcher.callHapB3HaplotypeMatches();
    }

    // check if any haplotype(s) show up on both strands in any diplotype
    Set<String> homozygous = new HashSet<>();
    for (DiplotypeMatch dm : matches) {
      Map<String, Integer> haps = new HashMap<>();
      for (String h : dm.getHaplotype1().getHaplotypeNames()) {
        if (!dm.getHaplotype1().getHaplotype().isPartial() || !h.equals(TextConstants.REFERENCE)) {
          haps.compute(h, (k, v) -> v == null ? 1 : v + 1);
        }
      }
      if (dm.getHaplotype2() != null) {
        for (String h : dm.getHaplotype2().getHaplotypeNames()) {
          if (!dm.getHaplotype2().getHaplotype().isPartial() || !h.equals(TextConstants.REFERENCE)) {
            haps.compute(h, (k, v) -> v == null ? 1 : v + 1);
          }
        }
      }
      for (String k : haps.keySet()) {
        if (haps.get(k) > 1) {
          homozygous.add(k);
        }
      }
    }

    // Reference/Reference matches should have been handled in the effectively phased section.
    // If we get to this point, it's a combination that might include Reference.
    if (dpydHapB3Matcher != null && dpydHapB3Matcher.isHapB3Present()) {
      // ...unless we have HapB3, which means it cannot be Reference/Reference
      homozygous.remove(TextConstants.REFERENCE);
    }

    SortedSet<HaplotypeMatch> hapMatches = comboData.comparePermutations();
    // If there are more than 2 haplotype matches, strip out Reference because we prioritize non-Reference if possible.
    // With 2 or fewer haplotype matches, cannot have Reference if there's a partial.
    int numMatches = hapMatches.size();
    if (dpydHapB3Matcher != null) {
      numMatches += dpydHapB3Matcher.getNumHapB3Called();
    }
    if (numMatches > 2) {
      hapMatches = hapMatches.stream()
          .filter(m -> !m.getName().equals(TextConstants.REFERENCE))
          .collect(Collectors.toCollection(TreeSet::new));
    }

    List<HaplotypeMatch> finalHaps = new ArrayList<>();
    for (HaplotypeMatch hm : hapMatches) {
      finalHaps.add(hm);
      if (homozygous.contains(hm.getName())) {
        finalHaps.add(hm);
        homozygous.remove(hm.getName());
      }
    }
    if (dpydHapB3Matcher != null) {
      finalHaps.addAll(dpydHapB3Matcher.buildHapB3HaplotypeMatches());
    }
    if (!homozygous.isEmpty()) {
      throw new IllegalStateException("Combination matching found " + homozygous + " but haplotype matching didn't");
    }
    return finalHaps;
  }


  /**
   * Initializes data required to call a diplotype.
   *
   * @param alleleMap map of {@link SampleAllele}s from VCF
   * @param assumeReference true if missing alleles in {@link NamedAllele}s should be treated as reference.
   * This determines how to interpret PharmCAT named allele definitions, not sample data.
   */
  private MatchData initializeCallData(String sampleId, SortedMap<String, SampleAllele> alleleMap, String gene,
      boolean assumeReference, boolean findCombinations) {

    SortedSet<NamedAllele> alleles = m_definitionReader.getHaplotypes(gene);
    VariantLocus[] allPositions = m_definitionReader.getPositions(gene);
    DefinitionExemption exemption = m_definitionReader.getExemption(gene);

    SortedSet<VariantLocus> extraPositions = null;
    if (exemption != null) {
      extraPositions = exemption.getExtraPositions();
    }

    // grab SampleAlleles for all positions related to the current gene
    MatchData data = new MatchData(sampleId, gene, alleleMap, allPositions, extraPositions, exemption);
    if (data.getNumSampleAlleles() == 0) {
      return data;
    }

    // handle missing positions (if any)
    data.marshallHaplotypes(gene, alleles, findCombinations);

    if (assumeReference) {
      // fill in blanks in named alleles based on the reference named allele
      data.defaultMissingAllelesToReference();
    }

    data.generateSamplePermutations();
    return data;
  }


  /**
   * Removes HapB3 NamedAlleles from {@link MatchData}.
   */
  private MatchData initializeDpydCallData(String sampleId, SortedMap<String, SampleAllele> alleleMap,
      boolean assumeReference, boolean findCombinations) {

    String gene = "DPYD";
    // remove HapB3 and HapB3Intron
    SortedSet<NamedAllele> alleles = m_definitionReader.getHaplotypes(gene).stream()
        .filter(a -> !a.getName().equals(DpydHapB3Matcher.HAPB3_ALLELE) &&
            !a.getName().equals(DpydHapB3Matcher.HAPB3_INTRONIC_ALLELE))
        .collect(Collectors.toCollection(TreeSet::new));
    VariantLocus[] allPositions = m_definitionReader.getPositions(gene);
    DefinitionExemption exemption = m_definitionReader.getExemption(gene);

    SortedSet<VariantLocus> extraPositions = null;
    if (exemption != null) {
      extraPositions = exemption.getExtraPositions();
    }

    // grab SampleAlleles for all positions related to the current gene
    MatchData data = new MatchData(sampleId, gene, alleleMap, allPositions, extraPositions, exemption);
    if (data.getNumSampleAlleles() == 0) {
      return data;
    }

    // handle missing positions (if any)
    data.marshallHaplotypes(gene, alleles, findCombinations);

    if (assumeReference) {
      data.defaultMissingAllelesToReference();
    }

    data.generateSamplePermutations();
    return data;
  }
}
