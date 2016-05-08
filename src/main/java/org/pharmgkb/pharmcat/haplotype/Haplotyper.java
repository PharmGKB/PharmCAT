package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Path;
import java.util.ArrayList;
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
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
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
  protected static ImmutableSet<String> calculateLocationsOfInterest(DefinitionReader definitionReader) {

    ImmutableSet.Builder<String> setBuilder = ImmutableSet.builder();
    for (String gene : definitionReader.getGenes()) {
      VariantLocus[] variants = definitionReader.getPositions(gene);
      String chromosome = definitionReader.getDefinitionFile(gene).getChromosome();
      // map variant to chr:position
      for (VariantLocus v : variants) {
        String chrPos = chromosome + ":" + v.getPosition();
        setBuilder.add(chrPos);
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
      report.gene(gene, callDiplotypes(alleles, gene), alleles.values());
    }
    return report;
  }


  /**
   * Calls the possible diplotypes for a single gene.
   *
   * @param alleleMap map of {@link SampleAllele}s from VCF
   */
  protected List<DiplotypeMatch> callDiplotypes(SortedMap<String, SampleAllele> alleleMap, String gene) {

    Dataset data = new Dataset();
    // grab SampleAlleles for all positions related to current gene
    data.marshallSampleData(alleleMap, m_definitionReader.getDefinitionFile(gene).getChromosome(),
        m_definitionReader.getPositions(gene));
    // handle missing positions (if any)
    data.marshallHaplotypes(m_definitionReader.getPositions(gene), m_definitionReader.getHaplotypes(gene));

    if (m_assumeReferenceInDefinitions) {
      data.defaultMissingAllelesToReference();
    }

    // get all permutations of sample at positions of interest
    Set<String> permutations = CombinationUtil.generatePermutations(
        data.geneSampleMap.values().stream()
            .sorted()
            .collect(Collectors.toList())
    );

    // find matched pairs
    List<DiplotypeMatch> pairs = new DiplotypeMatcher(data.geneSampleMap, permutations, data.haplotypes)
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


  static class Dataset {
    SortedMap<Integer, SampleAllele> geneSampleMap = new TreeMap<>();
    Set<VariantLocus> missingPositions = new HashSet<>();
    VariantLocus[] positions;
    List<NamedAllele> haplotypes;


    /**
     * Organizes the {@link SampleAllele} data related to the {@code gene} of interest.
     *
     * @param alleleMap map of {@link SampleAllele}s from VCF
     */
    void marshallSampleData(@Nonnull SortedMap<String, SampleAllele> alleleMap, @Nonnull String chromosome,
        @Nonnull VariantLocus[] positions) {

      for (VariantLocus variant : positions) {
        String chrPos = chromosome + ":" + variant.getPosition();
        SampleAllele allele = alleleMap.get(chrPos);
        if (allele == null) {
          missingPositions.add(variant);
          sf_logger.info("Sample has no allele for {}", chrPos);
          continue;
        }
        allele = allele.forVariant(variant);
        geneSampleMap.put(variant.getPosition(), allele);
      }
    }

    /**
     * Organizes the {@link NamedAllele} data for analysis.
     * This will also reorganize haplotypes to deal with samples that have missing alleles.
     */
    void marshallHaplotypes(VariantLocus[] allPositions, List<NamedAllele> allHaplotypes) {

      if (missingPositions.isEmpty()) {
        positions = allPositions;
        haplotypes = allHaplotypes;

      } else {
        // handle missing positions by duplicating haplotype and eliminating missing positions
        VariantLocus[] availablePositions = new VariantLocus[allPositions.length - missingPositions.size()];
        for (int x = 0, y = 0; x < allPositions.length; x += 1) {
          if (!missingPositions.contains(allPositions[x])) {
            availablePositions[y] = allPositions[x];
            y += 1;
          }
        }
        positions = availablePositions;
        haplotypes = new ArrayList<>();
        for (NamedAllele hap : allHaplotypes) {
          // get alleles for positions we have data on
          String[] availableAlleles = new String[positions.length];
          for (int x = 0; x < positions.length; x += 1) {
            availableAlleles[x] = hap.getAllele(positions[x]);
            if (availableAlleles[x] != null) {
              sf_logger.info("{} cannot be matched due to missing positions", hap.getName());
              break;
            }
          }
          NamedAllele newHap = new NamedAllele(hap.getId(), hap.getName(), availableAlleles);
          newHap.setFunction(hap.getFunction());
          newHap.setPopFreqMap(hap.getPopFreqMap());
          newHap.finalize(positions);
          if (newHap.getNumValidAlleles() > 0) {
            haplotypes.add(newHap);
          }
        }
      }
    }


    /**
     * Assumes that missing alleles in {@link NamedAllele}s should be the reference.
     */
    void defaultMissingAllelesToReference() {

      List<NamedAllele> updatedHaplotypes = new ArrayList<>();
      NamedAllele referenceHaplotype = null;
      for (NamedAllele hap : haplotypes) {
        if (referenceHaplotype == null) {
          referenceHaplotype = hap;
          updatedHaplotypes.add(hap);
          continue;
        }

        String[] refAlleles = referenceHaplotype.getAlleles();
        String[] curAlleles = hap.getAlleles();
        Preconditions.checkState(refAlleles.length == curAlleles.length);

        String[] newAlleles = new String[refAlleles.length];
        for (int x = 0; x < refAlleles.length; x += 1) {
          if (curAlleles[x] == null) {
            newAlleles[x] = refAlleles[x];
          } else {
            newAlleles[x] = curAlleles[x];
          }
        }

        NamedAllele fixedHap = new NamedAllele(hap.getId(), hap.getName(), newAlleles);
        fixedHap.setFunction(hap.getFunction());
        fixedHap.setPopFreqMap(hap.getPopFreqMap());
        fixedHap.finalize(positions);
        updatedHaplotypes.add(fixedHap);
      }

      haplotypes = updatedHaplotypes;
    }
  }
}
