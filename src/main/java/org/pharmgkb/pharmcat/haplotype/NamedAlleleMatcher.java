package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.concurrent.ThreadSafe;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableSet;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.annotation.AnnotationReader;
import org.pharmgkb.pharmcat.annotation.model.RsidAnnotation;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Result;


/**
 * This is the main entry point for matching {@link NamedAllele}s.
 *
 * @author Mark Woon
 */
@ThreadSafe
public class NamedAlleleMatcher {
  public static final String VERSION = "1.0.0";
  private DefinitionReader m_definitionReader;
  private AnnotationReader m_annotationReader;
  private ImmutableSet<String> m_locationsOfInterest;
  private boolean m_assumeReferenceInDefinitions;
  private boolean m_topCandidateOnly;


  /**
   * Default constructor.
   * This will only call the top candidate(s) and assume reference.
   */
  public NamedAlleleMatcher(@Nonnull DefinitionReader definitionReader, @Nonnull AnnotationReader annotationReader) {
    this(definitionReader, annotationReader, true, false);
  }

  /**
   * Constructor.
   *
   * @param topCandidateOnly true if only top candidate(s) should be called, false to call all possible candidates
   * @param assumeReference true if missing alleles in definitions should be treated as reference, false otherwise
   */
  public NamedAlleleMatcher(@Nonnull DefinitionReader definitionReader, @Nonnull AnnotationReader annotationReader,
      boolean assumeReference, boolean topCandidateOnly) {

    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_annotationReader = annotationReader;
    m_locationsOfInterest = calculateLocationsOfInterest(m_definitionReader);
    m_assumeReferenceInDefinitions = assumeReference;
    m_topCandidateOnly = topCandidateOnly;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("d", "definition-dir", "directory of allele definition files", true, "d")
          .addOption("a", "annotations-dir", "directory of annotation files", false, "a")
          .addOption("vcf", "vcf-in", "VCF file", true, "vcf")
          .addOption("json", "json-out", "file to save results to (in JSON format)", false, "json")
          .addOption("html", "html-out", "file to save results to (in HTML format)", false, "html");

      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path definitionDir = cliHelper.getValidDirectory("d", false);
      Path annotationDir = definitionDir;
      if (cliHelper.hasOption("a")) {
        annotationDir = cliHelper.getValidDirectory("a", false);
      }
      Path vcfFile = cliHelper.getValidFile("vcf", false);

      DefinitionReader definitionReader = new DefinitionReader();
      definitionReader.read(definitionDir);
      if (definitionReader.getGenes().size() == 0) {
        System.out.println("Did not find any allele definitions at " + definitionDir);
        System.exit(1);
      }

      AnnotationReader annotationReader = new AnnotationReader();
      annotationReader.read(annotationDir);

      NamedAlleleMatcher namedAlleleMatcher = new NamedAlleleMatcher(definitionReader, annotationReader);
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
  private static ImmutableSet<String> calculateLocationsOfInterest(DefinitionReader definitionReader) {

    ImmutableSet.Builder<String> setBuilder = ImmutableSet.builder();
    for (String gene : definitionReader.getGenes()) {
      Arrays.stream(definitionReader.getPositions(gene))
          .map(VariantLocus::getVcfChrPosition)
          .forEach(setBuilder::add);
    }
    return setBuilder.build();
  }


  /**
   * Calls diplotypes for the given VCF file for all genes for which a definition exists.
   */
  public Result call(@Nonnull Path vcfFile) throws IOException {

    VcfReader vcfReader = buildVcfReader(vcfFile);
    SortedMap<String, SampleAllele> alleles = vcfReader.getAlleleMap();
    ResultBuilder resultBuilder = new ResultBuilder(m_definitionReader)
        .forFile(vcfFile);
    // call haplotypes
    for (String gene : m_definitionReader.getGenes()) {
      MatchData data = initializeCallData(alleles, gene);

      List<DiplotypeMatch> matches = null;
      if (data.getNumSampleAlleles() > 0) {
        matches = callDiplotypes(data);
      }

      Set<SampleAllele> annotatedAlleles = new HashSet<>();
      Collection<RsidAnnotation> rsidAnnotations = m_annotationReader.getRsidAnnotations(gene);
      if (rsidAnnotations.size() > 0) {
        for (RsidAnnotation ann : rsidAnnotations) {
          VariantLocus varLoc = m_definitionReader.getDefinitionFile(ann.getGene()).getVariantByRsid(ann.getRsid());
          if (varLoc != null) {
            SampleAllele sa = alleles.get(varLoc.getVcfChrPosition());
            if (sa != null) {
              annotatedAlleles.add(sa);
            }
          }
        }
      }

      resultBuilder.gene(gene, data, matches, annotatedAlleles);
    }
    return resultBuilder.build();
  }


  /**
   * Initializes data required to call a diplotype.
   *
   * @param alleleMap map of {@link SampleAllele}s from VCF
   */
  private @Nonnull MatchData initializeCallData(SortedMap<String, SampleAllele> alleleMap, String gene) {

    // grab SampleAlleles for all positions related to current gene
    MatchData data = new MatchData(alleleMap, m_definitionReader.getPositions(gene));
    if (data.getNumSampleAlleles() == 0) {
      return data;
    }
    // handle missing positions (if any)
    data.marshallHaplotypes(m_definitionReader.getHaplotypes(gene));

    if (m_assumeReferenceInDefinitions) {
      data.defaultMissingAllelesToReference();
    }

    data.generateSamplePermutations();
    return data;
  }


  /**
   * Calls the possible diplotypes for a single gene.
   *
   */
  protected List<DiplotypeMatch> callDiplotypes(MatchData data) {

    // find matched pairs
    List<DiplotypeMatch> pairs = new DiplotypeMatcher(data)
        .compute();
    if (m_topCandidateOnly) {
      if (pairs.size() > 1) {
        int topScore = pairs.get(0).getScore();
        pairs = pairs.stream()
            .filter(dm -> dm.getScore() == topScore)
            .collect(Collectors.toList());
      }
    }
    return pairs;
  }
}
