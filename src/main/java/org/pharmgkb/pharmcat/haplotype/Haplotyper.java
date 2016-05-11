package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.SortedMap;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.concurrent.ThreadSafe;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableSet;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;


/**
 * This is the main entry point for calling haplotypes.
 *
 * @author Mark Woon
 */
@ThreadSafe
public class Haplotyper {
  private DefinitionReader m_definitionReader;
  private VcfReader m_vcfReader;
  private boolean m_assumeReferenceInDefinitions;
  private boolean m_topCandidateOnly;


  /**
   * Default constructor.
   * This will only call the top candidate(s) and assume reference.
   */
  public Haplotyper(@Nonnull DefinitionReader definitionReader) {
    this(definitionReader, true, true);
  }

  /**
   * Constructor.
   *
   * @param topCandidateOnly true if only top candidate(s) should be called, false to call all possible candidates
   * @param assumeReference true if missing alleles in definitions should be treated as reference, false otherwise
   */
  public Haplotyper(@Nonnull DefinitionReader definitionReader, boolean assumeReference, boolean topCandidateOnly) {

    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_vcfReader = new VcfReader(calculateLocationsOfInterest(m_definitionReader));
    m_assumeReferenceInDefinitions = assumeReference;
    m_topCandidateOnly = topCandidateOnly;
  }


  VcfReader getVcfReader() {
    return m_vcfReader;
  }


  /**
   * Collects all locations of interest (i.e. positions necessary to make a haplotype call).
   *
   * @return a set of {@code <chr:position>} Strings
   */
  private static ImmutableSet<String> calculateLocationsOfInterest(DefinitionReader definitionReader) {

    ImmutableSet.Builder<String> setBuilder = ImmutableSet.builder();
    for (String gene : definitionReader.getGenes()) {
      VariantLocus[] variants = definitionReader.getPositions(gene);
      String chromosome = definitionReader.getDefinitionFile(gene).getChromosome();
      // map variant to chr:position
      for (VariantLocus v : variants) {
        setBuilder.add(chromosome + ":" + v.getVcfPosition());
      }
    }
    return setBuilder.build();
  }


  /**
   * Calls diplotypes for the given VCF file for all genes for which a definition exists.
   */
  public Report call(@Nonnull Path vcfFile) throws IOException {

    SortedMap<String, SampleAllele> alleles = m_vcfReader.read(vcfFile);
    Report report = new Report(m_definitionReader)
        .forFile(vcfFile);
    // call haplotypes
    for (String gene : m_definitionReader.getGenes()) {
      MatchData data = initializeCallData(alleles, gene);
      report.gene(gene, callDiplotypes(data), data);
    }
    return report;
  }


  /**
   * Initializes data required to call a diplotype.
   *
   * @param alleleMap map of {@link SampleAllele}s from VCF
   */
  MatchData initializeCallData(SortedMap<String, SampleAllele> alleleMap, String gene) {

    MatchData data = new MatchData();
    // grab SampleAlleles for all positions related to current gene
    data.marshallSampleData(alleleMap, m_definitionReader.getDefinitionFile(gene).getChromosome(),
        m_definitionReader.getPositions(gene));
    if (data.geneSampleMap.size() == 0) {
      throw new IllegalStateException("No alleles in sample for " + gene);
    }
    // handle missing positions (if any)
    data.marshallHaplotypes(m_definitionReader.getPositions(gene), m_definitionReader.getHaplotypes(gene));

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
