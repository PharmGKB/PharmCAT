package org.pharmgkb.pharmcat.haplotype;

import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This is the data used to compute a {@link DiplotypeMatch} for a specific gene.
 *
 * @author Mark Woon
 */
public class MatchData {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private SortedMap<Integer, SampleAllele> m_sampleMap = new TreeMap<>();
  private VariantLocus[] m_positions;
  @Expose
  @SerializedName("missingPositions")
  private Set<VariantLocus> m_missingPositions = new HashSet<>();
  private List<NamedAllele> m_haplotypes;
  private Set<String> m_permutations;


  /**
   * Constructor.
   * Organizes the {@link SampleAllele} data related for the gene of interest.
   *
   * @param alleleMap map of chr:positions to {@link SampleAllele}s from VCF
   * @param allPositions all {@link VariantLocus} positions of interest for the gene
   */
  public MatchData(@Nonnull SortedMap<String, SampleAllele> alleleMap, @Nonnull String chromosome,
      @Nonnull VariantLocus[] allPositions) {

    List<VariantLocus> positions = new ArrayList<>();
    for (VariantLocus variant : allPositions) {
      String chrPos = chromosome + ":" + variant.getVcfPosition();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        m_missingPositions.add(variant);
        sf_logger.info("Sample has no allele for {}", chrPos);
        continue;
      }
      positions.add(variant);
      allele = allele.forVariant(variant);
      m_sampleMap.put(variant.getVcfPosition(), allele);
    }
    m_positions = positions.toArray(new VariantLocus[0]);
  }


  /**
   * Organizes the {@link NamedAllele} data for analysis.
   * This will also reorganize haplotypes to deal with samples that have missing alleles.
   */
  void marshallHaplotypes(List<NamedAllele> allHaplotypes) {

    if (m_missingPositions.isEmpty()) {
      m_haplotypes = allHaplotypes;

    } else {
      // handle missing positions by duplicating haplotype and eliminating missing positions
      m_haplotypes = new ArrayList<>();
      for (NamedAllele hap : allHaplotypes) {
        // get alleles for positions we have data on
        String[] availableAlleles = new String[m_positions.length];
        for (int x = 0; x < m_positions.length; x += 1) {
          availableAlleles[x] = hap.getAllele(m_positions[x]);
        }

        SortedSet<VariantLocus> missingPositions = m_missingPositions.stream()
            .filter(l -> hap.getAllele(l) != null)
            .collect(Collectors.toCollection(TreeSet::new));

        NamedAllele newHap = new NamedAllele(hap.getId(), hap.getName(), availableAlleles, missingPositions);
        newHap.setFunction(hap.getFunction());
        newHap.setPopFreqMap(hap.getPopFreqMap());
        newHap.initialize(m_positions);
        if (newHap.getScore() > 0) {
          m_haplotypes.add(newHap);
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
    for (NamedAllele hap : m_haplotypes) {
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

      NamedAllele fixedHap = new NamedAllele(hap.getId(), hap.getName(), newAlleles, hap.getMissingPositions());
      fixedHap.setFunction(hap.getFunction());
      fixedHap.setPopFreqMap(hap.getPopFreqMap());
      fixedHap.initialize(m_positions, hap.getScore());
      updatedHaplotypes.add(fixedHap);
    }

    m_haplotypes = updatedHaplotypes;
  }


  public int getNumSampleAlleles() {
    return m_sampleMap.size();
  }

  public @Nonnull SampleAllele getSampleAllele(int position) {
    SampleAllele sampleAllele = m_sampleMap.get(position);
    if (sampleAllele == null) {
      throw new IllegalArgumentException("No sample allele for position " + position);
    }
    return sampleAllele;
  }

  /**
   * Gets all permutations of sample alleles at positions of interest.
   */
  public @Nonnull Set<String> getPermutations() {
    if (m_permutations == null) {
      throw new IllegalStateException("Not initialized - call generateSamplePermutations()");
    }
    return m_permutations;
  }

  /**
   * Generate all permutations of sample alleles at positions of interest.
   */
  void generateSamplePermutations() {

    m_permutations = CombinationUtil.generatePermutations(
        m_sampleMap.values().stream()
            .sorted()
            .collect(Collectors.toList())
    );
  }


  /**
   * Gets the positions available for calling the haplotypes for the gene.
   */
  public @Nonnull VariantLocus[] getPositions() {
    return m_positions;
  }

  /**
   * Gets the positions that are missing from the sample VCF that would have been helpful for calling the haplotypes for
   * the gene.
   */
  public @Nonnull Set<VariantLocus> getMissingPositions() {
    return m_missingPositions;
  }

  /**
   * Gets the callable haplotypes for the gene based on the available positions.
   */
  public @Nonnull List<NamedAllele> getHaplotypes() {
    if (m_haplotypes == null) {
      if (m_sampleMap.size() == 0) {
        return Collections.emptyList();
      }
      throw new IllegalStateException("Not initialized - call marshallHaplotypes()");
    }
    return m_haplotypes;
  }
}
