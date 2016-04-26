package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.concurrent.ThreadSafe;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableSet;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This is the main entry point for calling haplotypes.
 *
 * @author Mark Woon
 */
@ThreadSafe
public class Haplotyper {
  private static Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private DefinitionReader m_definitionReader;
  private VcfReader m_vcfReader;
  private boolean m_topCandidateOnly;


  /**
   * Default constructor.
   * Will only call top candidate(s) by default.
   */
  public Haplotyper(@Nonnull DefinitionReader definitionReader) {
    this(definitionReader, true);
  }

  /**
   * Constructor.
   *
   * @param topCandidateOnly true if only top candidate(s) should be called, false to call all possible candidates
   */
  public Haplotyper(@Nonnull DefinitionReader definitionReader, boolean topCandidateOnly) {

    Preconditions.checkNotNull(definitionReader);
    m_definitionReader = definitionReader;
    m_vcfReader = new VcfReader(calculateLocationsOfInterest(m_definitionReader));
    m_topCandidateOnly = topCandidateOnly;
  }


  /**
   * Collects all locations of interest (i.e. positions necessary to make a haplotype call).
   *
   * @return a set of {@code <chr:position>} Strings
   */
  protected static ImmutableSet<String> calculateLocationsOfInterest(DefinitionReader definitionReader) {

    ImmutableSet.Builder<String> setBuilder = ImmutableSet.builder();
    for (String gene : definitionReader.getHaplotypes().keySet()) {
      List<Variant> variants = definitionReader.getHaplotypePositions().get(gene);
      // map variant to chr:position
      for (Variant v : variants) {
        String chrPos = v.getChromosome() + ":" + v.getPosition();
        setBuilder.add(chrPos);
      }
      // calculate permutations of all haplotypes
      for (Haplotype hap : definitionReader.getHaplotypes().get(gene)) {
        hap.calculatePermutations(variants);
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
    for (String gene : m_definitionReader.getHaplotypePositions().keySet()) {
      report.gene(gene, callDiplotypes(alleles, gene), alleles.values());
    }
    return report;
  }


  /**
   * Calls the possible diplotypes for a single gene.
   */
  protected List<DiplotypeMatch> callDiplotypes(SortedMap<String, SampleAllele> alleleMap, String gene) {

    // grab SampleAlleles for all positions related to current gene
    SortedMap<Integer, SampleAllele> geneSampleMap = new TreeMap<>();
    Set<String> missingPositions = new HashSet<>();
    for (Variant variant : m_definitionReader.getHaplotypePositions().get(gene)) {
      String chrPos = variant.getChromosome() + ":" + variant.getPosition();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        missingPositions.add(chrPos);
        sf_logger.info("Sample has no allele for {} (ref is {})", chrPos, variant.getREF());
        continue;
      }
      geneSampleMap.put(variant.getPosition(), allele);
    }

    // get all permutations of sample at positions of interest
    Set<String> permutations = CombinationUtil.generatePermutations(
        geneSampleMap.values().stream()
            .sorted()
            .collect(Collectors.toList())
    );


    // handle missing positions of interest in sample
    List<Haplotype> haplotypes;
    if (missingPositions.isEmpty()) {
      haplotypes = m_definitionReader.getHaplotypes().get(gene);
    } else {
      // TODO(markwoon): can probably just clone haplotype here and recompute permutations
      throw new UnsupportedOperationException("Sample is missing alleles for " + missingPositions);
    }

    // find matched pairs
    List<DiplotypeMatch> pairs = new DiplotypeMatcher(geneSampleMap, permutations, haplotypes).compute();
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
