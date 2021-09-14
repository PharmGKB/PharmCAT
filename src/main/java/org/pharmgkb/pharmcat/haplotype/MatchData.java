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
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This is the data used to compute a {@link DiplotypeMatch} for a specific gene.
 *
 * @author Mark Woon
 */
public class MatchData {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private final SortedMap<Integer, SampleAllele> m_sampleMap = new TreeMap<>();
  private final VariantLocus[] m_positions;
  @Expose
  @SerializedName("missingPositions")
  private final SortedSet<VariantLocus> m_missingPositions = new TreeSet<>();
  private final Set<VariantLocus> m_ignoredPositions = new HashSet<>();
  private final SortedSet<Variant> m_extraPositions = new TreeSet<>();
  @Expose
  @SerializedName("mismatchedAlleles")
  private final SortedSet<VariantLocus> m_mismatchedAlleles = new TreeSet<>();
  private SortedSet<NamedAllele> m_haplotypes;
  private Set<String> m_permutations;


  /**
   * Constructor.
   * Organizes the {@link SampleAllele} data related for the gene of interest.
   *
   * @param alleleMap map of chr:positions to {@link SampleAllele}s from VCF
   * @param allPositions all {@link VariantLocus} positions of interest for the gene
   * @param extraPositions extra positions to track sample alleles for
   * @param ignoredPositions ignored positions due to ignored named alleles
   */
  public MatchData(SortedMap<String, SampleAllele> alleleMap, VariantLocus[] allPositions,
      @Nullable SortedSet<VariantLocus> extraPositions, @Nullable SortedSet<VariantLocus> ignoredPositions) {

    if (ignoredPositions != null) {
      m_ignoredPositions.addAll(ignoredPositions);
    }

    List<VariantLocus> positions = new ArrayList<>();
    for (VariantLocus variant : allPositions) {
      String chrPos = variant.getVcfChrPosition();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        m_missingPositions.add(variant);
        sf_logger.debug("Sample has no allele for {}", chrPos);
        continue;
      }
      if (m_ignoredPositions.contains(variant)) {
        continue;
      }
      positions.add(variant);
      allele = allele.forVariant(variant);
      m_sampleMap.put(variant.getVcfPosition(), allele);
    }
    m_positions = positions.toArray(new VariantLocus[0]);
    if (extraPositions != null) {
      for (VariantLocus vl : extraPositions) {
        SampleAllele allele = alleleMap.get(vl.getVcfChrPosition());
        if (allele != null) {
          m_extraPositions.add(new Variant(vl, allele));
        } else {
          m_extraPositions.add(new Variant(vl.getPosition(), vl.getRsid(), null, vl.getVcfPosition(), null));
        }
      }
    }
  }


  void checkAlleles(DefinitionFile definitionFile) {

    for (VariantLocus variantLocus : m_positions) {
      SampleAllele sampleAllele = m_sampleMap.get(variantLocus.getVcfPosition());
      if (sampleAllele != null) {
        Set<String> alleles = definitionFile.getVariantAlleles(variantLocus);
        if (!alleles.contains(sampleAllele.getAllele1()) ||
            (sampleAllele.getAllele2() != null && !alleles.contains(sampleAllele.getAllele2()))) {
          m_mismatchedAlleles.add(variantLocus);
        }
      }
    }
  }


  /**
   * Organizes the {@link NamedAllele} data for analysis.
   * This will also reorganize haplotypes to deal with samples that have missing alleles.
   */
  void marshallHaplotypes(SortedSet<NamedAllele> allHaplotypes) {

    if (m_missingPositions.isEmpty() && m_ignoredPositions.isEmpty()) {
      m_haplotypes = allHaplotypes;

    } else {
      // handle missing positions by duplicating haplotype and eliminating missing positions
      m_haplotypes = new TreeSet<>();
      for (NamedAllele hap : allHaplotypes) {
        // get alleles for positions we have data on
        String[] availableAlleles = new String[m_positions.length];
        for (int x = 0; x < m_positions.length; x += 1) {
          availableAlleles[x] = hap.getAllele(m_positions[x]);
        }

        SortedSet<VariantLocus> missingPositions = m_missingPositions.stream()
            .filter(l -> hap.getAllele(l) != null)
            .collect(Collectors.toCollection(TreeSet::new));

        NamedAllele newHap = new NamedAllele(hap.getId(), hap.getName(), availableAlleles, missingPositions, hap.isReference());
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

    SortedSet<NamedAllele> updatedHaplotypes = new TreeSet<>();
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

      NamedAllele fixedHap = new NamedAllele(hap.getId(), hap.getName(), newAlleles, hap.getMissingPositions(), hap.isReference());
      fixedHap.setPopFreqMap(hap.getPopFreqMap());
      fixedHap.initialize(m_positions, hap.getScore());
      updatedHaplotypes.add(fixedHap);
    }

    m_haplotypes = updatedHaplotypes;
  }


  public int getNumSampleAlleles() {
    return m_sampleMap.size();
  }

  public SampleAllele getSampleAllele(int position) {
    SampleAllele sampleAllele = m_sampleMap.get(position);
    if (sampleAllele == null) {
      throw new IllegalArgumentException("No sample allele for position " + position);
    }
    return sampleAllele;
  }

  /**
   * Gets all permutations of sample alleles at positions of interest.
   */
  public Set<String> getPermutations() {
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
  public VariantLocus[] getPositions() {
    return m_positions;
  }

  /**
   * Gets the positions that are missing from the sample VCF that would have been helpful for calling the haplotypes for
   * the gene.
   */
  public SortedSet<VariantLocus> getMissingPositions() {
    return m_missingPositions;
  }

  /**
   * Gets the positions that are mismatched from any allele defined for the given gene (i.e. at given position, no
   * alleles in VCF match what we expect to see).
   *
   * @return a Set of {@link VariantLocus} objects with mismatched alleles
   */
  public Set<VariantLocus> getMismatchedPositions() {
    return m_mismatchedAlleles;
  }

  /**
   * Gets the extra positions specified in {@link DefinitionExemption#getExtraPositions()}.
   */
  public SortedSet<Variant> getExtraPositions() {
    return m_extraPositions;
  }

  /**
   * Gets the callable haplotypes for the gene based on the available positions.
   */
  public SortedSet<NamedAllele> getHaplotypes() {
    if (m_haplotypes == null) {
      if (m_sampleMap.size() == 0) {
        return Collections.emptySortedSet();
      }
      throw new IllegalStateException("Not initialized - call marshallHaplotypes()");
    }
    return m_haplotypes;
  }
}
