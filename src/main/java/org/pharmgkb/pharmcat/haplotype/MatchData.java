package org.pharmgkb.pharmcat.haplotype;

import java.lang.invoke.MethodHandles;
import java.util.*;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
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
  private final String m_sampleId;
  private final String m_gene;
  private final SortedMap<Long, SampleAllele> m_sampleMap = new TreeMap<>();
  private final boolean m_isHaploid;
  /** Positions at which data is available for sample. */
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
  @Expose
  @SerializedName("phased")
  private boolean m_isPhased;
  @Expose
  @SerializedName("homozygous")
  private boolean m_isHomozygous;
  @Expose
  @SerializedName("effectivelyPhased")
  private boolean m_isEffectivelyPhased;
  @Expose
  @SerializedName("hasNovelAllele")
  private boolean m_hasNovelAllele;
  private final Map<String, Map<Object, Object>> m_sequenceAlleleCache = new HashMap<>();


  /**
   * Constructor.
   * Organizes the {@link SampleAllele} data related for the gene of interest.
   *
   * @param alleleMap map of chr:positions to {@link SampleAllele}s from VCF
   * @param allPositions all {@link VariantLocus} positions of interest for the gene
   * @param extraPositions extra positions to track sample alleles for
   * @param ignoredPositions ignored positions due to ignored named alleles
   */
  public MatchData(String sampleId, String gene, SortedMap<String, SampleAllele> alleleMap, VariantLocus[] allPositions,
      @Nullable SortedSet<VariantLocus> extraPositions, @Nullable SortedSet<VariantLocus> ignoredPositions) {
    m_sampleId = sampleId;
    m_gene = gene;
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
      if (allele.getNovelAlleles().size() > 0) {
        m_hasNovelAllele = true;
      }
      positions.add(variant);
      m_sampleMap.put(variant.getPosition(), allele);
    }
    m_positions = positions.toArray(new VariantLocus[0]);
    if (extraPositions != null) {
      for (VariantLocus vl : extraPositions) {
        SampleAllele allele = alleleMap.get(vl.getVcfChrPosition());
        if (allele != null) {
          m_extraPositions.add(new Variant(vl, allele));
        } else {
          m_extraPositions.add(new Variant(vl.getPosition(), vl.getRsid(), null, null));
        }
      }
    }
    m_isHaploid = m_sampleMap.values().stream().allMatch(sa -> sa.getAllele2() == null);
    m_isPhased = m_sampleMap.values().stream().allMatch(SampleAllele::isPhased);
    m_isHomozygous = m_isHaploid ||
        m_sampleMap.values().stream().allMatch(sa -> sa.getAllele1().equals(sa.getAllele2()));
  }


  void checkAlleles() {

    for (VariantLocus variantLocus : m_positions) {
      SampleAllele sampleAllele = m_sampleMap.get(variantLocus.getPosition());
      if (sampleAllele != null) {
        if (!variantLocus.hasVcfAllele(sampleAllele.getAllele1()) ||
            (sampleAllele.getAllele2() != null && !variantLocus.hasVcfAllele(sampleAllele.getAllele2()))) {
          m_mismatchedAlleles.add(variantLocus);
        }
      }
    }
  }


  /**
   * Organizes the {@link NamedAllele} data for analysis.
   * This will also reorganize haplotypes to deal with samples that have missing alleles.
   */
  void marshallHaplotypes(String gene, SortedSet<NamedAllele> allHaplotypes, boolean findCombinations) {

    if (m_missingPositions.isEmpty() && m_ignoredPositions.isEmpty()) {
      if (findCombinations) {
        m_haplotypes = new TreeSet<>();
        for (NamedAllele hap : allHaplotypes) {
          if (isIgnorableCombination(gene, hap)) {
            continue;
          }
          m_haplotypes.add(hap);
        }
      } else {
        m_haplotypes = allHaplotypes;
      }

    } else {
      // handle missing positions by duplicating haplotype and eliminating missing positions
      m_haplotypes = new TreeSet<>();
      for (NamedAllele hap : allHaplotypes) {
        if (findCombinations) {
          if (isIgnorableCombination(gene, hap)) {
            continue;
          }
        }
        // get alleles for positions we have data on
        String[] availableAlleles = new String[m_positions.length];
        String[] cpicAlleles = new String[m_positions.length];
        for (int x = 0; x < m_positions.length; x += 1) {
          availableAlleles[x] = hap.getAllele(m_positions[x]);
          cpicAlleles[x] = hap.getCpicAllele(m_positions[x]);
        }

        SortedSet<VariantLocus> missingPositions = m_missingPositions.stream()
            .filter(l -> hap.getAllele(l) != null)
            .collect(Collectors.toCollection(TreeSet::new));

        NamedAllele newHap = new NamedAllele(hap.getId(), hap.getName(), availableAlleles, cpicAlleles,
            missingPositions, hap.isReference());
        newHap.setPopFreqMap(hap.getPopFreqMap());
        newHap.initialize(m_positions);
        if (newHap.getScore() > 0) {
          m_haplotypes.add(newHap);
        }
      }
    }
  }

  private boolean isIgnorableCombination(String gene, NamedAllele hap) {
    if (gene.equalsIgnoreCase("UGT1A1")) {
      return hap.getName().contains("+");
    }
    return false;
  }


  /**
   * Assumes that missing alleles in {@link NamedAllele}s should be the reference.
   */
  void defaultMissingAllelesToReference() {

    SortedSet<NamedAllele> updatedHaplotypes = new TreeSet<>();
    Optional<NamedAllele> refHapOpt = m_haplotypes.stream().filter(NamedAllele::isReference).findAny();
    if (refHapOpt.isEmpty()) {
      throw new IllegalStateException(m_gene + " does not have a reference");
    }
    NamedAllele referenceHaplotype = refHapOpt.get();
    int numAlleles = referenceHaplotype.getAlleles().length;
    for (NamedAllele hap : m_haplotypes) {
      if (referenceHaplotype == hap) {
        updatedHaplotypes.add(hap);
        continue;
      }

      String[] curAlleles = hap.getAlleles();
      Preconditions.checkState(numAlleles == curAlleles.length);

      String[] newAlleles = new String[numAlleles];
      String[] cpicAlleles = new String[numAlleles];
      for (int x = 0; x < numAlleles; x += 1) {
        if (curAlleles[x] == null) {
          // ref allele can be null if position is missing
          String refAllele = referenceHaplotype.getAllele(x);
          if (refAllele != null && Iupac.isWobble(refAllele)) {
            newAlleles[x] = m_positions[x].getRef();
          } else {
            newAlleles[x] = refAllele;
          }
          cpicAlleles[x] = referenceHaplotype.getCpicAlleles()[x];
        } else {
          newAlleles[x] = curAlleles[x];
          cpicAlleles[x] = hap.getCpicAlleles()[x];
        }
      }

      NamedAllele fixedHap = new NamedAllele(hap.getId(), hap.getName(), newAlleles, cpicAlleles,
          hap.getMissingPositions(), hap.isReference());
      fixedHap.setPopFreqMap(hap.getPopFreqMap());
      fixedHap.initialize(m_positions, hap.getScore());
      updatedHaplotypes.add(fixedHap);
    }

    m_haplotypes = updatedHaplotypes;
  }


  public int getNumSampleAlleles() {
    return m_sampleMap.size();
  }

  public SampleAllele getSampleAllele(long position) {
    SampleAllele sampleAllele = m_sampleMap.get(position);
    if (sampleAllele == null) {
      throw new IllegalArgumentException("No sample allele for position " + position);
    }
    return sampleAllele;
  }

  public boolean isHaploid() {
    return m_isHaploid;
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
            .toList()
    );
    m_isEffectivelyPhased = m_permutations.size() <= 2;
  }

  /**
   * Gets whether data is phased (i.e. is phased at all positions).
   */
  public boolean isPhased() {
    return m_isPhased;
  }

  public boolean isHomozygous() {
    return m_isHomozygous;
  }

  /**
   * Gets whether data is "effectively phased" (i.e. actually phased or unphased but homozygous at all positions).
   */
  public boolean isEffectivelyPhased() {
    return m_isEffectivelyPhased;
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


  /**
   * Utility method to cache allele lookups in sequences.
   */
  public String getAllele(String sequence, int idx) {
    Map<Object, Object> seqMap = m_sequenceAlleleCache.computeIfAbsent(sequence, s -> {
      Map<Object, Object> m = new HashMap<>();
      m.put("all", s.split(";"));
      return m;
    });

    String allele = (String)seqMap.get(idx);
    if (allele == null) {
      allele = ((String[])seqMap.get("all"))[idx].split(":")[1];
      seqMap.put(idx, allele);
    }
    return allele;
  }


  public boolean hasNovelAllele() {
    return m_hasNovelAllele;
  }


  /**
   * Compares a sample's allele permutations to haplotype definitions and return matches.
   */
  protected SortedSet<HaplotypeMatch> comparePermutations() {
    Set<HaplotypeMatch> haplotypeMatches = getHaplotypes().stream()
        .map(HaplotypeMatch::new)
        .collect(Collectors.toSet());
    for (String p : getPermutations()) {
      for (HaplotypeMatch hm : haplotypeMatches) {
        hm.match(p);
      }
    }
    return haplotypeMatches.stream()
        .filter(h -> !h.getSequences().isEmpty())
        .collect(Collectors.toCollection(TreeSet::new));
  }


  @Override
  public String toString() {
    return m_gene + " match data for " + m_sampleId;
  }
}
