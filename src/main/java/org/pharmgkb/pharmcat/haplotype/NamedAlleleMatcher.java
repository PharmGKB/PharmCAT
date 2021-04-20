package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
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
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * This is the main entry point for matching {@link NamedAllele}s.
 *
 * @author Mark Woon
 */
public class NamedAlleleMatcher {
  public static final String VERSION = "1.0.0";
  private final DefinitionReader m_definitionReader;
  private final ImmutableMap<String, VariantLocus> m_locationsOfInterest;
  private final boolean m_assumeReferenceInDefinitions;
  private final boolean m_topCandidateOnly;
  private boolean m_printWarnings;


  /**
   * Default constructor.
   * This will only call the top candidate(s) and assume reference.
   */
  public NamedAlleleMatcher(DefinitionReader definitionReader) {
    this(definitionReader, true, false);
  }

  /**
   * Constructor.
   *
   * @param topCandidateOnly true if only top candidate(s) should be called, false to call all possible candidates
   * @param assumeReference true if missing alleles in definitions should be treated as reference, false otherwise
   */
  public NamedAlleleMatcher(DefinitionReader definitionReader, boolean assumeReference,
      boolean topCandidateOnly) {

    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_locationsOfInterest = calculateLocationsOfInterest(m_definitionReader);
    m_assumeReferenceInDefinitions = assumeReference;
    m_topCandidateOnly = topCandidateOnly;
  }


  public NamedAlleleMatcher printWarnings() {
    m_printWarnings = true;
    return this;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("vcf", "vcf-in", "VCF file", true, "vcf")
          .addOption("json", "json-out", "file to save results to (in JSON format)", false, "json")
          .addOption("html", "html-out", "file to save results to (in HTML format)", false, "html")
          .addOption("d", "definition-dir", "directory of allele definition files", false, "d")
          .addOption("a", "all-results", "return all possible results, not just top hits")
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
      NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(definitionReader, true, topCandidateOnly)
          .printWarnings();
      Result result = namedAlleleMatcher.call(vcfFile);

      ResultSerializer resultSerializer = new ResultSerializer();
      if (cliHelper.hasOption("json")) {
        resultSerializer.toJson(result, cliHelper.getPath("json"));
      }
      if (cliHelper.hasOption("html")) {
        resultSerializer.toHtml(result, cliHelper.getPath("html"));
      }

    } catch (Exception ex) {
      ex.printStackTrace();
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
  private static ImmutableMap<String, VariantLocus> calculateLocationsOfInterest(DefinitionReader definitionReader) {

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


  /**
   * Calls diplotypes for the given VCF file for all genes for which a definition exists.
   */
  public Result call(Path vcfFile) throws IOException {

    VcfReader vcfReader = buildVcfReader(vcfFile);
    SortedMap<String, SampleAllele> alleles = vcfReader.getAlleleMap();
    ResultBuilder resultBuilder = new ResultBuilder(m_definitionReader)
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
      DefinitionExemption exemption = m_definitionReader.getExemption(gene);
      MatchData data = initializeCallData(alleles, gene);
      List<DiplotypeMatch> matches = null;
      if (data.getNumSampleAlleles() > 0) {
        boolean topCandidateOnly = m_topCandidateOnly;
        if (exemption != null && exemption.isAllHits() != null) {
          topCandidateOnly = !exemption.isAllHits();
        }
        matches = callDiplotypes(data, topCandidateOnly);
      }

      resultBuilder.gene(gene, data, matches);
    }
    return resultBuilder.build();
  }


  /**
   * Initializes data required to call a diplotype.
   *
   * @param alleleMap map of {@link SampleAllele}s from VCF
   */
  private MatchData initializeCallData(SortedMap<String, SampleAllele> alleleMap, String gene) {

    DefinitionExemption exemption = m_definitionReader.getExemption(gene);
    SortedSet<VariantLocus> extraPositions = null;
    List<NamedAllele> alleles = m_definitionReader.getHaplotypes(gene);
    VariantLocus[] allPositions = m_definitionReader.getPositions(gene);
    SortedSet<VariantLocus> unusedPositions = null;
    if (exemption != null) {
      extraPositions = exemption.getExtraPositions();
      unusedPositions = findUnusedPositions(exemption, allPositions, alleles);
    }

    // grab SampleAlleles for all positions related to current gene
    MatchData data = new MatchData(alleleMap, allPositions, extraPositions, unusedPositions);
    data.checkAlleles(m_definitionReader.getDefinitionFile(gene));
    if (data.getNumSampleAlleles() == 0) {
      return data;
    }

    if (exemption != null) {
      alleles = alleles.stream()
          .filter(a -> !exemption.shouldIgnoreAllele(a.getName()))
          .collect(Collectors.toList());
    }
    // handle missing positions (if any)
    data.marshallHaplotypes(alleles);

    boolean assumeReference = m_assumeReferenceInDefinitions;
    if (exemption != null && exemption.isAssumeReference() != null) {
      assumeReference = exemption.isAssumeReference();
    }
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
      List<NamedAllele> namedAlleles) {

    SortedSet<VariantLocus> unusedPositions = new TreeSet<>();
    if (exemption.getIgnoredAlleles().isEmpty()) {
      return unusedPositions;
    }

    List<NamedAllele> variantNamedAlleles = namedAlleles.subList(1, namedAlleles.size() - 1);
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


  /**
   * Calls the possible diplotypes for a single gene.
   *
   */
  protected List<DiplotypeMatch> callDiplotypes(MatchData data, boolean topCandidateOnly) {

    // find matched pairs
    List<DiplotypeMatch> pairs = new DiplotypeMatcher(data)
        .compute();
    if (topCandidateOnly && pairs.size() > 1) {
      int topScore = pairs.get(0).getScore();
      pairs = pairs.stream()
          .filter(dm -> dm.getScore() == topScore)
          .collect(Collectors.toList());
    }
    return pairs;
  }
}
