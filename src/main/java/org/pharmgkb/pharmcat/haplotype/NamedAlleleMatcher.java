package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
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
import com.google.common.collect.Sets;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.BaseConfig;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.DataManager;


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
      if (gene.equals("DPYD")) {
        callDpyd(vcfReader.getSampleId(), alleleMap, resultBuilder);
      } else if (gene.equals("RYR1")) {
        callLowestFunctionGene(vcfReader.getSampleId(), "RYR1", alleleMap, resultBuilder);
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

    MatchData comboData = null;
    // try for diplotypes if effectively phased
    if (origData.isEffectivelyPhased()) {
      // first look for exact matches (use topCandidateOnly = false because looking for exact match)
      SortedSet<DiplotypeMatch> diplotypeMatches = new DiplotypeMatcher(m_env, origData)
          .compute(false, false);
      if (diplotypeMatches.size() == 1) {
        resultBuilder.diplotypes(gene, origData, diplotypeMatches);
        return;
      }
      // try combinations
      comboData = initializeCallData(sampleId, alleleMap, gene, false, true);
      diplotypeMatches = callPhasedLowestFunctionGeneWithCombination(comboData, origData.isHomozygous());
      if (!diplotypeMatches.isEmpty()) {
        resultBuilder.diplotypes(gene, origData, diplotypeMatches);
        return;
      }
    }

    if (comboData == null) {
      comboData = initializeCallData(sampleId, alleleMap, gene, false, true);
    }
    List<HaplotypeMatch> hapMatches = callHaplotypesForLowestFunctionGene(comboData, null);
    resultBuilder.haplotypes(gene, origData, hapMatches);
  }


  private SortedSet<DiplotypeMatch> callPhasedLowestFunctionGeneWithCombination(MatchData comboData,
      boolean isHomozygous) {

    SortedSet<DiplotypeMatch> diplotypeMatches = new DiplotypeMatcher(m_env, comboData)
        .compute(true, false, false);
    if (!diplotypeMatches.isEmpty()) {
      DiplotypeMatch[] matches = diplotypeMatches.toArray(new DiplotypeMatch[0]);
      for (int x = 0; x < matches.length; x += 1) {
        matches = removeSubCombos(matches, x);
      }
      SortedSet<DiplotypeMatch> finalMatches = new TreeSet<>(Arrays.asList(matches));
      if (matches.length > 1) {
        if (isHomozygous) {
          finalMatches = finalMatches.stream()
              .filter(m -> m.getHaplotype2() != null &&
                  m.getHaplotype1().getName().equals(m.getHaplotype2().getName()))
              .collect(Collectors.toCollection(TreeSet::new));
        }
      }
      if (finalMatches.size() > 1) {
        // this should never happen
        throw new IllegalStateException("Lowest function gene cannot have more than 1 diplotype");
      }
      return finalMatches;
    }
    return diplotypeMatches;
  }


  private List<HaplotypeMatch> callHaplotypesForLowestFunctionGene(MatchData comboData,
      @Nullable DpydHapB3Matcher dpydHapB3Matcher) {

    SortedSet<HaplotypeMatch> hapMatches = comboData.comparePermutations();
    // have to compute diplotypes so that we can check for homozygous and partials
    SortedSet<DiplotypeMatch> matches = new DiplotypeMatcher(m_env, comboData)
        .compute(true, false, false);
    Set<String> homozygous = new HashSet<>();
    int numPartials = 0;
    for (DiplotypeMatch dm : matches) {
      if (dm.getHaplotype1().getHaplotype().isPartial() ||
          (dm.getHaplotype2() != null && dm.getHaplotype2().getHaplotype().isPartial())) {
        numPartials += 1;
      }

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

    // If there are more than 2 haplotype matches, strip out Reference because we prioritize non-Reference if possible.
    // With 2 or fewer haplotype matches, cannot have Reference if there's a partial.
    int numMatches = hapMatches.size();
    if (dpydHapB3Matcher != null) {
      numMatches += dpydHapB3Matcher.getNumHapB3Called();
    }
    if (numMatches > 2 || numPartials > 0) {
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
    if (!homozygous.isEmpty()) {
      throw new IllegalStateException("Combination matching found " + homozygous + " but haplotype matching didn't");
    }
    return finalHaps;
  }


  /**
   * For DPYD, only attempt an exact diplotype match if phased.
   * If there is no exact match, we only look for potential haplotypes.
   * This tries to match all permutations to any potential haplotype (won't assume reference).
   */
  private void callDpyd(String sampleId, SortedMap<String, SampleAllele> alleleMap, ResultBuilder resultBuilder) {
    final String gene = "DPYD";

    MatchData origData = initializeCallData(sampleId, alleleMap, gene, true, false);
    if (origData.getNumSampleAlleles() == 0) {
      resultBuilder.noCall(gene, origData);
      return;
    }
    if (!hasRequiredPositions(origData, resultBuilder)) {
      return;
    }

    DpydHapB3Matcher dpydHapB3Matcher = new DpydHapB3Matcher(m_env, alleleMap, origData.isEffectivelyPhased());
    MatchData workingData;
    if (dpydHapB3Matcher.isMissingHapB3Positions()) {
      workingData = origData;
    } else {
      workingData = initializeDpydCallData(sampleId, alleleMap, true, false);
    }

    MatchData comboData = null;
    // try for diplotypes if effectively phased
    if (origData.isEffectivelyPhased()) {
     // first look for exact matches (use topCandidateOnly = false because looking for exact match)
      SortedSet<DiplotypeMatch> diplotypeMatches = new DiplotypeMatcher(m_env, workingData)
          .compute(false, false);
      if (diplotypeMatches.size() == 1) {
        if (!dpydHapB3Matcher.isMissingHapB3Positions()) {
          MatchData mergerData = initializeCallData(sampleId, alleleMap, gene, false, false);
          SortedSet<DiplotypeMatch> mergedMatches = dpydHapB3Matcher.mergePhasedDiplotypeMatch(mergerData, diplotypeMatches);
          resultBuilder.diplotypes(gene, mergerData, mergedMatches, dpydHapB3Matcher.getWarnings());
        } else {
          resultBuilder.diplotypes(gene, workingData, diplotypeMatches);
        }
        return;
      }
      // try combinations
      comboData = initializeDpydCallData(sampleId, alleleMap, false, true);
      diplotypeMatches = callPhasedLowestFunctionGeneWithCombination(comboData, workingData.isHomozygous());
      if (!diplotypeMatches.isEmpty()) {
        if (!dpydHapB3Matcher.isMissingHapB3Positions()) {
          MatchData mergerData = initializeCallData(sampleId, alleleMap, gene, false, true);
          SortedSet<DiplotypeMatch> mergedMatches = dpydHapB3Matcher.mergePhasedDiplotypeMatch(mergerData, diplotypeMatches);
          resultBuilder.diplotypes(gene, mergerData, mergedMatches, dpydHapB3Matcher.getWarnings());
        } else {
          resultBuilder.diplotypes(gene, comboData, diplotypeMatches);
        }
        return;
      }
    }

    if (comboData == null) {
      comboData = initializeDpydCallData(sampleId, alleleMap, false, true);
    }
    List<HaplotypeMatch> hapMatches = callHaplotypesForLowestFunctionGene(comboData, dpydHapB3Matcher);
    if (dpydHapB3Matcher.isHapB3Present()) {
      hapMatches.addAll(dpydHapB3Matcher.buildHapB3HaplotypeMatches(origData));
    }
    resultBuilder.haplotypes(gene, origData, hapMatches, dpydHapB3Matcher.getWarnings());
  }

  /**
   * Removes diplotypes that are partial combinations.
   */
  private DiplotypeMatch[] removeSubCombos(DiplotypeMatch[] diplotypeMatches, int idx) {

    List<DiplotypeMatch> cleaned = new ArrayList<>(Arrays.asList(diplotypeMatches).subList(0, idx));
    DiplotypeMatch diplotypeMatch = diplotypeMatches[idx];
    Set<String> hap1 = getNames(diplotypeMatch.getHaplotype1());
    Set<String> hap2 = getNames(diplotypeMatch.getHaplotype2());
    cleaned.add(diplotypeMatch);
    for (int x = idx + 1; x < diplotypeMatches.length; x += 1) {
      DiplotypeMatch dm = diplotypeMatches[x];
      Set<String> h1 = getNames(dm.getHaplotype1());
      Set<String> h2 = getNames(dm.getHaplotype2());
      if ((!hap1.containsAll(h1) || !hap2.containsAll(h2)) &&
          (!hap1.containsAll(h2) || !hap2.containsAll(h1))) {
        cleaned.add(dm);
      }
    }
    return cleaned.toArray(new DiplotypeMatch[0]);
  }

  private Set<String> getNames(@Nullable BaseMatch bm) {
    if (bm == null) {
      return Collections.emptySet();
    }
    if (bm instanceof CombinationMatch cm) {
      return cm.getComponentHaplotypes().stream()
          .map(NamedAllele::getName)
          .collect(Collectors.toSet());
    } else {
      return Sets.newHashSet(bm.getName());
    }
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
    MatchData data = new MatchData(sampleId, gene, alleleMap, allPositions, extraPositions, null, exemption);
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
    SortedSet<VariantLocus> unusedPositions = new TreeSet<>();
    // add HapB3 positions to ignore
    for (VariantLocus vl : allPositions) {
      if (DpydHapB3Matcher.isHapB3Rsid(vl.getRsid())) {
        unusedPositions.add(vl);
      }
    }

    // grab SampleAlleles for all positions related to the current gene
    MatchData data = new MatchData(sampleId, gene, alleleMap, allPositions, extraPositions, unusedPositions, exemption);
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
