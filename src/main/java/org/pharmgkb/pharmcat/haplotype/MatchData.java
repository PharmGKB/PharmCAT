package org.pharmgkb.pharmcat.haplotype;

import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This is the data used to compute a {@link DiplotypeMatch}.
 *
 * @author Mark Woon
 */
public class MatchData {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private SortedMap<Integer, SampleAllele> m_sampleMap = new TreeMap<>();
  private Set<VariantLocus> m_missingPositions = new HashSet<>();
  private VariantLocus[] m_positions;
  private List<NamedAllele> m_haplotypes;
  private Set<String> m_permutations;


  /**
   * Constructor.
   * Organizes the {@link SampleAllele} data related for a single gene of interest.
   *
   * @param alleleMap map of chr:positions to {@link SampleAllele}s from VCF
   */
  public MatchData(@Nonnull SortedMap<String, SampleAllele> alleleMap, @Nonnull String chromosome,
      @Nonnull VariantLocus[] positions) {

    for (VariantLocus variant : positions) {
      String chrPos = chromosome + ":" + variant.getVcfPosition();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        m_missingPositions.add(variant);
        sf_logger.info("Sample has no allele for {}", chrPos);
        continue;
      }
      allele = allele.forVariant(variant);
      m_sampleMap.put(variant.getVcfPosition(), allele);
    }
  }


  /**
   * Organizes the {@link NamedAllele} data for analysis.
   * This will also reorganize haplotypes to deal with samples that have missing alleles.
   */
  void marshallHaplotypes(VariantLocus[] allPositions, List<NamedAllele> allHaplotypes) {

    if (m_missingPositions.isEmpty()) {
      m_positions = allPositions;
      m_haplotypes = allHaplotypes;

    } else {
      // handle missing positions by duplicating haplotype and eliminating missing positions
      VariantLocus[] availablePositions = new VariantLocus[allPositions.length - m_missingPositions.size()];
      for (int x = 0, y = 0; x < allPositions.length; x += 1) {
        if (!m_missingPositions.contains(allPositions[x])) {
          availablePositions[y] = allPositions[x];
          y += 1;
        }
      }
      m_positions = availablePositions;
      m_haplotypes = new ArrayList<>();
      for (NamedAllele hap : allHaplotypes) {
        // get alleles for positions we have data on
        String[] availableAlleles = new String[m_positions.length];
        for (int x = 0; x < m_positions.length; x += 1) {
          availableAlleles[x] = hap.getAllele(m_positions[x]);
        }
        NamedAllele newHap = new NamedAllele(hap.getId(), hap.getName(), availableAlleles);
        newHap.setFunction(hap.getFunction());
        newHap.setPopFreqMap(hap.getPopFreqMap());
        newHap.finalize(m_positions);
        if (newHap.getNumValidAlleles() > 0) {
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

      NamedAllele fixedHap = new NamedAllele(hap.getId(), hap.getName(), newAlleles);
      fixedHap.setFunction(hap.getFunction());
      fixedHap.setPopFreqMap(hap.getPopFreqMap());
      fixedHap.finalize(m_positions);
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

  public @Nonnull Set<String> getPermutations() {
    if (m_positions == null) {
      throw new IllegalStateException("Not initialized - call generateSamplePermutations()");
    }
    return m_permutations;
  }

  /**
   * Generate all permutations of sample at positions of interest.
   */
  void generateSamplePermutations() {

    m_permutations = CombinationUtil.generatePermutations(
        m_sampleMap.values().stream()
            .sorted()
            .collect(Collectors.toList())
    );
  }


  public @Nonnull VariantLocus[] getPositions() {
    if (m_positions == null) {
      throw new IllegalStateException("Not initialized - call marshallHaplotypes()");
    }
    return m_positions;
  }

  public @Nonnull Set<VariantLocus> getMissingPositions() {
    return m_missingPositions;
  }

  public @Nonnull List<NamedAllele> getHaplotypes() {
    if (m_haplotypes == null) {
      throw new IllegalStateException("Not initialized - call marshallHaplotypes()");
    }
    return m_haplotypes;
  }
}
