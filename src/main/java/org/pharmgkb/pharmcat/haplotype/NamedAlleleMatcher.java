package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Sets;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * This is the main entry point for matching {@link NamedAllele}s.
 *
 * @author Mark Woon
 */
public class NamedAlleleMatcher {
  public static final String VERSION = "2.0.0";
  private final DefinitionReader m_definitionReader;
  private final ImmutableMap<String, VariantLocus> m_locationsOfInterest;
  private final boolean m_findCombinations;
  private final boolean m_topCandidateOnly;
  private final boolean m_callCyp2d6;
  private boolean m_printWarnings;


  /**
   * Default constructor.
   * This will only call the top candidate(s).  Will not find combinations or call CYP2D6.
   */
  public NamedAlleleMatcher(DefinitionReader definitionReader) {
    this(definitionReader, false, false, false);
  }

  /**
   * Constructor.
   *
   * @param topCandidateOnly true if only top candidate(s) should be called, false to call all possible candidates
   * @param callCyp2d6 true if CYP2D6 should be called
   */
  public NamedAlleleMatcher(DefinitionReader definitionReader, boolean findCombinations, boolean topCandidateOnly,
      boolean callCyp2d6) {

    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_locationsOfInterest = calculateLocationsOfInterest(m_definitionReader);
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

      DefinitionReader definitionReader = new DefinitionReader();
      definitionReader.read(definitionDir);
      if (definitionReader.getGenes().size() == 0) {
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
          new NamedAlleleMatcher(definitionReader, findCombinations, topCandidateOnly, callCyp2d6)
              .printWarnings();
      Result result = namedAlleleMatcher.call(vcfFile);

      Path jsonFile = CliUtils.getOutputFile(cliHelper, vcfFile, "json", ".match.json");
      ResultSerializer resultSerializer = new ResultSerializer();
      resultSerializer.toJson(result, jsonFile);
      System.out.println("Saved JSON results to " + jsonFile);
      if (cliHelper.hasOption("html")) {
        Path htmlFile = CliUtils.getOutputFile(cliHelper, vcfFile, "html", ".match.html");
        resultSerializer.toHtml(result, htmlFile);
        System.out.println("Saved HTML results to " + htmlFile);
      }

    } catch (Exception ex) {
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
   * Builds a new VCF reader for the given file.
   */
  VcfReader buildVcfReader(Path vcfFile) throws IOException {
    return new VcfReader(m_locationsOfInterest, vcfFile);
  }


  /**
   * Collects all locations of interest (i.e. positions necessary to make a haplotype call).
   *
   * @return a set of {@code <chr:position>} Strings
   */
  static ImmutableMap<String, VariantLocus> calculateLocationsOfInterest(DefinitionReader definitionReader) {

    Set<String> data = new HashSet<>();
    ImmutableMap.Builder<String, VariantLocus> mapBuilder = ImmutableMap.builder();
    for (String gene : definitionReader.getGenes()) {
      Arrays.stream(definitionReader.getPositions(gene))
          .forEach(v -> {
            String vcp = v.getVcfChrPosition();
            data.add(vcp);
            mapBuilder.put(vcp, v);
          });
      DefinitionExemption exemption = definitionReader.getExemption(gene);
      if (exemption != null) {
        exemption.getExtraPositions()
            .forEach(v -> {
              String vcp = v.getVcfChrPosition();
              if (!data.contains(vcp)) {
                data.add(vcp);
                mapBuilder.put(vcp, v);
              }
            });
      }
    }
    return mapBuilder.build();
  }


  private boolean getTopCandidateOnly(String gene) {
    DefinitionExemption exemption = m_definitionReader.getExemption(gene);
    if (exemption != null && exemption.isAllHits() != null) {
      //noinspection ConstantConditions
      return !exemption.isAllHits();
    }
    return m_topCandidateOnly;
  }


  /**
   * Calls diplotypes for the given VCF file for all genes for which a definition exists.
   */
  public Result call(Path vcfFile) throws IOException {

    VcfReader vcfReader = buildVcfReader(vcfFile);
    SortedMap<String, SampleAllele> alleleMap = vcfReader.getAlleleMap();
    ResultBuilder resultBuilder = new ResultBuilder(m_definitionReader, m_topCandidateOnly, m_findCombinations, m_callCyp2d6)
        .forFile(vcfFile, vcfReader.getWarnings().asMap());
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
        callDpyd(alleleMap, resultBuilder);
      } else {
        callAssumingReference(alleleMap, gene, resultBuilder);
      }
    }
    return resultBuilder.build();
  }

  private void callAssumingReference(SortedMap<String, SampleAllele> alleleMap, String gene,
      ResultBuilder resultBuilder) {

    MatchData data = initializeCallData(alleleMap, gene, true, false);
    if (data.getNumSampleAlleles() == 0) {
      resultBuilder.gene(gene, data);
      return;
    }

    List<DiplotypeMatch> matches = new DiplotypeMatcher(data)
        .compute(false, getTopCandidateOnly(gene));
    if (matches.isEmpty()) {
      if (!m_findCombinations) {
        resultBuilder.gene(gene, data, matches);
      } else {
        callCombination(alleleMap, gene, resultBuilder);
      }
      return;
    }
    resultBuilder.gene(gene, data, matches);
  }

  private void callCombination(SortedMap<String, SampleAllele> alleleMap, String gene, ResultBuilder resultBuilder) {

    MatchData data = initializeCallData(alleleMap, gene, false, true);
    if (data.getNumSampleAlleles() == 0) {
      resultBuilder.gene(gene, data);
      return;
    }

    List<DiplotypeMatch> matches = new DiplotypeMatcher(data)
        .compute(true, getTopCandidateOnly(gene));
    resultBuilder.gene(gene, data, matches);
  }


  /**
   * For DPYD, only attempt exact diplotype match if phased.
   * If there is no exact match, we only look for potential haplotypes.
   * This tries to match all permutations to any potential haplotype (won't assume reference).
   */
  private void callDpyd(SortedMap<String, SampleAllele> alleleMap, ResultBuilder resultBuilder) {
    final String gene = "DPYD";

    MatchData refData = initializeCallData(alleleMap, gene, true, false);
    if (refData.getNumSampleAlleles() == 0) {
      resultBuilder.gene(gene, refData);
      return;
    }

    MatchData comboData = null;
    // try for diplotypes if effectively phased
    if (refData.isEffectivelyPhased()) {
     // first look for exact matches (use topCandidateOnly = false because looking for exact match)
      List<DiplotypeMatch> diplotypeMatches = new DiplotypeMatcher(refData)
          .compute(false, false);
      if (diplotypeMatches.size() == 1) {
        resultBuilder.gene(gene, refData, diplotypeMatches);
        return;
      }
      // try combinations
      comboData = initializeCallData(alleleMap, gene, false, true);
      diplotypeMatches = new DiplotypeMatcher(comboData)
          .compute(true, false, false, true);
      if (!diplotypeMatches.isEmpty()) {
        DiplotypeMatch[] matches = diplotypeMatches.toArray(new DiplotypeMatch[0]);
        for (int x = 0; x < matches.length; x += 1) {
          matches = removeSubCombos(matches, x);
        }
        resultBuilder.gene(gene, comboData, Arrays.asList(matches));
        return;
      }
    }

    if (comboData == null) {
      comboData = initializeCallData(alleleMap, gene, false, true);
    }
    SortedSet<HaplotypeMatch> hapMatches = comboData.comparePermutations();
    if (hapMatches.size() > 2) {
      // cannot have reference if there is more than 1 variant
      hapMatches = hapMatches.stream()
          .filter(m -> !m.getName().equals("Reference"))
          .collect(Collectors.toCollection(TreeSet::new));
    }
    resultBuilder.gene(gene, refData, hapMatches);
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

  private Set<String> getNames(BaseMatch bm) {
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
   */
  private MatchData initializeCallData(SortedMap<String, SampleAllele> alleleMap, String gene,
      boolean assumeReference, boolean findCombinations) {

    DefinitionExemption exemption = m_definitionReader.getExemption(gene);
    SortedSet<VariantLocus> extraPositions = null;
    SortedSet<NamedAllele> alleles = m_definitionReader.getHaplotypes(gene);
    VariantLocus[] allPositions = m_definitionReader.getPositions(gene);
    SortedSet<VariantLocus> unusedPositions = null;
    if (exemption != null) {
      extraPositions = exemption.getExtraPositions();
      unusedPositions = findUnusedPositions(exemption, allPositions, alleles);
    }

    // grab SampleAlleles for all positions related to current gene
    MatchData data = new MatchData(alleleMap, allPositions, extraPositions, unusedPositions);
    data.checkAlleles();
    if (data.getNumSampleAlleles() == 0) {
      return data;
    }

    if (exemption != null) {
      alleles = alleles.stream()
          .filter(a -> !exemption.shouldIgnoreAllele(a.getName()))
          .collect(Collectors.toCollection(TreeSet::new));
    }
    // handle missing positions (if any)
    data.marshallHaplotypes(gene, alleles, findCombinations);

    if (assumeReference) {
      data.defaultMissingAllelesToReference();
    }

    data.generateSamplePermutations();
    return data;
  }


  /**
   * Find positions that are only used by ignored alleles (and therefore should be eliminated from consideration).
   */
  private SortedSet<VariantLocus> findUnusedPositions(DefinitionExemption exemption, VariantLocus[] allPositions,
      SortedSet<NamedAllele> namedAlleles) {

    SortedSet<VariantLocus> unusedPositions = new TreeSet<>();
    if (exemption.getIgnoredAlleles().isEmpty()) {
      return unusedPositions;
    }

    List<NamedAllele> allAlleles = new ArrayList<>(namedAlleles);
    List<NamedAllele> variantNamedAlleles = allAlleles.subList(1, namedAlleles.size() - 1);
    Set<VariantLocus> ignorablePositions = new HashSet<>();
    for (NamedAllele namedAllele : variantNamedAlleles) {
      if (exemption.shouldIgnoreAllele(namedAllele.getName())) {
        ignorablePositions.addAll(findIgnorablePositions(allPositions, namedAllele));
      }
    }

    for (VariantLocus vl : ignorablePositions) {
      boolean isUnused = true;
      for (NamedAllele namedAllele : variantNamedAlleles) {
        if (!exemption.shouldIgnoreAllele(namedAllele.getName())) {
          if (namedAllele.getAllele(vl) != null) {
            isUnused = false;
            break;
          }
        }
      }
      if (isUnused) {
        unusedPositions.add(vl);
      }
    }
    return unusedPositions;
  }

  /**
   * Find positions that are used by ignored alleles (and are therefore potentially ignoreable).
   */
  private Set<VariantLocus> findIgnorablePositions(VariantLocus[] allPositions, NamedAllele namedAllele)  {
    Set<VariantLocus> ignorablePositions = new HashSet<>();
    int x = 0;
    for (String allele : namedAllele.getAlleles()) {
      if (allele != null) {
        ignorablePositions.add(allPositions[x]);
      }
      x += 1;
    }
    return ignorablePositions;
  }
}
