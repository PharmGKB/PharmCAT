package org.pharmgkb.pharmcat.haplotype;

import java.lang.invoke.MethodHandles;
import java.util.*;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.Variant;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Mutable per-sample, per-gene working data used to compute {@link DiplotypeMatch}es.
 *
 * <p>The normal lifecycle is: construct, {@link #marshallHaplotypes(String, SortedSet, boolean)}, optionally
 * {@link #defaultMissingAllelesToReference()}, {@link #generateSamplePermutations()}, then
 * {@link #comparePermutations()}. Marshalling or reference defaulting can replace the callable haplotypes, so both
 * operations invalidate the lazily built candidate index.</p>
 *
 * @author Mark Woon
 */
public class MatchData {
  public static final Integer NULL_PHASE_SET = Integer.MIN_VALUE;
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private final String m_sampleId;
  private final String m_gene;
  private final SortedMap<Long, SampleAllele> m_sampleMap = new TreeMap<>();
  private final boolean m_isHaploid;
  /** Positions available for this sample, retained in definition/output order. */
  private final VariantLocus[] m_positions;
  /** Positions in the order used by internal sample permutations. */
  private final VariantLocus[] m_permutationPositions;
  /** Maps VCF position to the index of the matching VariantLocus in the m_positions array. */
  private final SortedMap<Long, Integer> m_vcfPositionIndex = new TreeMap<>();
  /** Maps VCF position to the index used by internal sample permutations. */
  private final Map<Long, Integer> m_permutationPositionIndex = new HashMap<>();
  @Expose
  @SerializedName("missingPositions")
  private final SortedSet<VariantLocus> m_missingPositions = new TreeSet<>();
  private final SortedSet<Long> m_missingVcfPositions = new TreeSet<>();
  private final SortedSet<Variant> m_extraPositions = new TreeSet<>();
  @Expose
  @SerializedName("positionsWithUndocumentedVariations")
  private final SortedSet<VariantLocus> m_positionsWithUndocumentedVariations = new TreeSet<>();
  @Expose
  @SerializedName("treatUndocumentedVariationsAsReference")
  private boolean m_treatUndocumentedVariationsAsReference;
  /** Sample-adjusted named alleles that remain callable at the available positions. */
  private @Nullable SortedSet<NamedAllele> m_haplotypes;
  /** Structured sample strands; their allele arrays follow m_permutationPositions order. */
  private @Nullable Set<SamplePermutation> m_permutations;
  /** Legacy/result representation of m_permutations, created only when requested. */
  private @Nullable Set<String> m_encodedPermutations;
  /** True when reference defaults are resolved lazily rather than copied into every named allele. */
  private boolean m_defaultMissingAllelesToReference;
  @Expose
  @SerializedName("phased")
  private final boolean m_isPhased;
  /** Map of phase set to positions. */
  @Expose
  @SerializedName("phaseSets")
  private final SortedMap<Integer, SortedSet<Long>> m_phaseSets = new TreeMap<>();
  /** Map of positions to phase set. */
  @Expose
  @SerializedName("posToPhaseSet")
  private final SortedMap<Long, Integer> m_positionToPhaseSet = new TreeMap<>();
  @Expose
  @SerializedName("homozygous")
  private boolean m_isHomozygous;
  @Expose
  @SerializedName("effectivelyPhased")
  private boolean m_isEffectivelyPhased;
  /** Parsed legacy sequences, used only when structured SamplePermutation metadata is unavailable. */
  private final Map<String, String[]> m_sequenceAlleleCache = new HashMap<>();
  /** Stable named-allele list whose indexes are the bit positions used by the candidate index. */
  private @Nullable List<NamedAllele> m_haplotypeIndex;
  /** Expected alleles aligned with m_haplotypeIndex and m_permutationPositions. */
  private @Nullable List<@Nullable String[]> m_haplotypeAlleles;
  /** Lazily materialized reference-defaulted objects aligned with m_haplotypeIndex for result output. */
  private @Nullable List<@Nullable NamedAllele> m_outputHaplotypes;
  /** Per permutation position: observed allele to compatible haplotype bits. Entries are populated on demand. */
  private @Nullable List<Map<String, BitSet>> m_candidateIndex;
  @Expose
  @SerializedName("missingRequiredPositions")
  private final List<String> m_missingRequiredPositions = new ArrayList<>();
  @Expose
  @SerializedName("missingAmp1Positions")
  private final List<String> m_missingAmp1Positions = new ArrayList<>();


  /**
   * Constructor.
   * Organizes the {@link SampleAllele} data related for the gene of interest.
   *
   * @param alleleMap map of chr:positions to {@link SampleAllele}s from VCF
   * @param allPositions all {@link VariantLocus} positions of interest for the gene
   * @param extraPositions extra positions to track sample alleles for
   */
  public MatchData(String sampleId, String gene, SortedMap<String, SampleAllele> alleleMap, VariantLocus[] allPositions,
      @Nullable SortedSet<VariantLocus> extraPositions, @Nullable DefinitionExemption exemption) {
    m_sampleId = sampleId;
    m_gene = gene;

    List<VariantLocus> positions = new ArrayList<>();
    boolean isPhased = true;
    for (VariantLocus variant : allPositions) {
      String chrPos = variant.getVcfChrPosition();
      SampleAllele allele = alleleMap.get(chrPos);
      if (allele == null) {
        m_missingPositions.add(variant);
        m_missingVcfPositions.add(variant.getPosition());
        sf_logger.debug("Sample has no allele for {}", chrPos);
        continue;
      }
      if (!allele.getUndocumentedVariations().isEmpty()) {
        m_positionsWithUndocumentedVariations.add(variant);
        if (allele.isTreatUndocumentedVariationsAsReference()) {
          m_treatUndocumentedVariationsAsReference = true;
        }
      }
      positions.add(variant);
      if (!allele.isPhased()) {
        isPhased = false;
      }
      if (allele.isPhased()) {
        m_phaseSets.computeIfAbsent(Objects.requireNonNullElse(allele.getPhaseSet(), NULL_PHASE_SET),
                ps -> new TreeSet<>())
            .add((long)allele.getPosition());
      }
      m_sampleMap.put(variant.getPosition(), allele);
    }
    m_positions = positions.toArray(new VariantLocus[0]);
    for (int x = 0; x < m_positions.length; x += 1) {
      m_vcfPositionIndex.put(m_positions[x].getPosition(), x);
    }
    m_permutationPositions = Arrays.stream(m_positions).sorted().toArray(VariantLocus[]::new);
    for (int x = 0; x < m_permutationPositions.length; x += 1) {
      m_permutationPositionIndex.put(m_permutationPositions[x].getPosition(), x);
    }
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
    m_isHaploid = areSampleAllelesHaploid(m_sampleMap.values());
    m_isPhased = isPhased;
    m_isHomozygous = m_isHaploid ||
        m_sampleMap.values().stream().allMatch(SampleAllele::isHomozygous);

    if (isUsingPhaseSets()) {
      for (Integer ps : m_phaseSets.keySet()) {
        for (Long pos : m_phaseSets.get(ps)) {
          m_positionToPhaseSet.put(pos, ps);
        }
      }
    }

    if (exemption != null && !m_missingPositions.isEmpty()) {
      if (exemption.hasRequiredPositions()) {
        for (VariantLocus missing : m_missingPositions) {
          if (exemption.isRequiredPosition(missing.getPosition())) {
            m_missingRequiredPositions.add(missing.getVcfChrPosition());
          }
        }
      }
      if (exemption.hasAmp1Positions()) {
        for (VariantLocus missing : m_missingPositions) {
          if (exemption.isAmp1Position(missing.getPosition())) {
            m_missingAmp1Positions.add(missing.getVcfChrPosition());
          }
        }
      }
    }
  }


  public String getSampleId() {
    return m_sampleId;
  }

  public String getGene() {
    return m_gene;
  }


  /**
   * Builds the callable {@link NamedAllele} set for the positions present in this sample.
   * When positions are missing, this creates sample-specific definitions with those positions removed and drops
   * definitions that no longer have a positive score.
   */
  void marshallHaplotypes(String gene, SortedSet<NamedAllele> allHaplotypes, boolean findCombinations) {

    m_defaultMissingAllelesToReference = false;
    clearCandidateIndex();

    if (m_missingPositions.isEmpty()) {
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
        @Nullable String[] availableAlleles = new String[m_positions.length];
        @Nullable String[] cpicAlleles = new String[m_positions.length];
        for (int x = 0; x < m_positions.length; x += 1) {
          availableAlleles[x] = hap.getAllele(m_positions[x]);
          cpicAlleles[x] = hap.getCpicAllele(m_positions[x]);
        }

        SortedSet<VariantLocus> missingPositions = m_missingPositions.stream()
            .filter(l -> hap.getAllele(l) != null)
            .collect(Collectors.toCollection(TreeSet::new));

        NamedAllele newHap = new NamedAllele(hap.getId(), hap.getName(), availableAlleles, cpicAlleles,
            missingPositions, hap.isReference());
        newHap.initialize(m_positions);
        if (newHap.getScore() > 0) {
          m_haplotypes.add(newHap);
        }
      }
    }
  }

  /**
   * Checks if any of the sample's alleles is partially missing.
   */
  public boolean hasPartialMissingAlleles() {
    return m_sampleMap.values().stream()
        .anyMatch(sa -> sa.getVcfCall().contains("."));
  }

  private boolean isIgnorableCombination(String gene, NamedAllele hap) {
    if (gene.equalsIgnoreCase("UGT1A1")) {
      return hap.getName().contains("+");
    }
    return false;
  }

  /**
   * Gets if this dataset is missing a required position.
   */
  public List<String> getMissingRequiredPositions() {
    return m_missingRequiredPositions;
  }

  /**
   * Gets if this dataset is missing a required position for AMP 1.
   */
  public List<String> getMissingAmp1Positions() {
    return m_missingAmp1Positions;
  }

  /**
   * Interprets blank cells in non-reference {@link NamedAllele} definitions as reference alleles.
   * With complete position data, defaults remain lazy and are applied by the candidate index. With missing position
   * data, the sample-specific definitions created by {@link #marshallHaplotypes(String, SortedSet, boolean)} are
   * replaced eagerly.
   */
  void defaultMissingAllelesToReference() {
    if (m_haplotypes == null) {
      throw new IllegalStateException("Not initialized - call marshallHaplotypes()");
    }

    NamedAllele referenceHaplotype = m_haplotypes.stream().filter(NamedAllele::isReference).findAny()
        .orElseThrow(() -> new IllegalStateException(m_gene + " does not have a reference"));
    if (m_missingPositions.isEmpty()) {
      m_defaultMissingAllelesToReference = true;
      clearCandidateIndex();
      return;
    }

    SortedSet<NamedAllele> updatedHaplotypes = new TreeSet<>();
    int numAlleles = referenceHaplotype.getAlleles().length;
    for (NamedAllele hap : m_haplotypes) {
      if (referenceHaplotype == hap) {
        updatedHaplotypes.add(hap);
        continue;
      }

      @Nullable String[] curAlleles = hap.getAlleles();
      Preconditions.checkState(numAlleles == curAlleles.length);

      @Nullable String[] newAlleles = new String[numAlleles];
      @Nullable String[] cpicAlleles = new String[numAlleles];
      for (int x = 0; x < numAlleles; x += 1) {
        if (curAlleles[x] == null) {
          // ref allele can be null if the position is missing
          String refAllele = referenceHaplotype.getAllele(x);
          if (Iupac.isWobble(refAllele)) {
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
      fixedHap.initialize(m_positions, hap.getScore());
      updatedHaplotypes.add(fixedHap);
    }

    m_haplotypes = updatedHaplotypes;
    m_defaultMissingAllelesToReference = false;
    clearCandidateIndex();
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
  Set<SamplePermutation> getPermutations() {
    if (m_permutations == null) {
      throw new IllegalStateException("Not initialized - call generateSamplePermutations()");
    }
    return m_permutations;
  }

  int getPermutationCount() {
    return getPermutations().size();
  }

  public Set<String> getPermutationStrings() {
    if (m_permutations == null) {
      throw new IllegalStateException("Not initialized - call generateSamplePermutations()");
    }
    if (m_encodedPermutations == null) {
      m_encodedPermutations = new HashSet<>();
      for (SamplePermutation permutation : m_permutations) {
        m_encodedPermutations.add(permutation.getSequence());
      }
    }
    return m_encodedPermutations;
  }


  /**
   * Generates all possible sample strands at positions of interest.
   * {@link CombinationUtil} receives sample alleles in numeric position order, which is also the order expected by
   * the structured matching path.
   */
  void generateSamplePermutations() {

    m_permutations = CombinationUtil.generatePermutationData(
        m_sampleMap.values().stream()
            .sorted()
            .toList()
    );
    m_encodedPermutations = null;
    m_isEffectivelyPhased = m_permutations.size() <= 2;
  }

  /**
   * Gets whether data is phased (i.e. is phased at all positions).
   */
  public boolean isPhased() {
    return m_isPhased;
  }

  /**
   * Gets whether data uses phase sets.
   */
  public boolean isUsingPhaseSets() {
    if (m_phaseSets.isEmpty()) {
      return false;
    }
    if (m_phaseSets.containsKey(NULL_PHASE_SET)) {
      return m_phaseSets.size() > 1;
    }
    return true;
  }

  /**
   * Gets a map of the positions for each phase set (i.e. phase set ID to positions).
   */
  public SortedMap<Integer, SortedSet<Long>> getPhaseSets() {
    return m_phaseSets;
  }

  /**
   * Gets the phase set ID for the specified {@code position}.
   */
  public @Nullable Integer getPhaseSet(long position) {
    return m_positionToPhaseSet.get(position);
  }


  public boolean isHomozygous() {
    return m_isHomozygous;
  }

  /**
   * Gets whether data is "effectively phased" (i.e., actually phased or unphased but homozygous at all positions).
   * More specifically, there is a maximum of 2 permutations of this sample's alleles.
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
   * Gets the positions that are missing from the sample VCF that would have been helpful for calling the haplotypes for
   * the gene.
   */
  public SortedSet<Long> getMissingVcfPositions() {
    return m_missingVcfPositions;
  }

  /**
   * Gets the positions that have variations that are not documented in the allele definition (i.e. any ALT alleles in
   * VCF that do not match what we expect to see).
   *
   * @return a Set of {@link VariantLocus} objects with undocumented variations
   */
  public Set<VariantLocus> getPositionsWithUndocumentedVariations() {
    return m_positionsWithUndocumentedVariations;
  }

  public boolean isTreatUndocumentedVariationsAsReference() {
    return m_treatUndocumentedVariationsAsReference;
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
      if (m_sampleMap.isEmpty()) {
        return Collections.emptySortedSet();
      }
      throw new IllegalStateException("Not initialized - call marshallHaplotypes()");
    }
    return m_haplotypes;
  }

  /**
   * Gets callable haplotypes with lazy reference defaults materialized for serialization and debug rendering.
   */
  SortedSet<NamedAllele> getHaplotypesForOutput() {
    SortedSet<NamedAllele> haplotypes = getHaplotypes();
    if (!m_defaultMissingAllelesToReference) {
      return haplotypes;
    }

    SortedSet<NamedAllele> outputHaplotypes = new TreeSet<>();
    NamedAllele referenceHaplotype = getReferenceHaplotype();
    for (NamedAllele haplotype : haplotypes) {
      if (haplotype.isReference()) {
        outputHaplotypes.add(haplotype);
      } else {
        outputHaplotypes.add(materializeDefaultedHaplotype(haplotype, referenceHaplotype));
      }
    }
    return outputHaplotypes;
  }


  /**
   * Utility method to cache allele lookups in sequences.
   */
  public String getAllele(String sequence, int idx) {
    return getSequenceAlleles(sequence)[idx];
  }

  public @Nullable String getAllele(String sequence, long vcfPosition) {
    Integer idx = m_vcfPositionIndex.get(vcfPosition);
    if (idx == null) {
      return null;
    }
    return getSequenceAlleles(sequence)[idx];
  }

  @Nullable String getAllele(SamplePermutation permutation, long vcfPosition) {
    Integer idx = m_permutationPositionIndex.get(vcfPosition);
    if (idx == null) {
      return null;
    }
    return permutation.getAllelesForMatching()[idx];
  }

  String[] getSequenceAlleles(String sequence) {
    return m_sequenceAlleleCache.computeIfAbsent(sequence, s -> {
      String[] alleles = new String[m_positions.length];
      for (String positionAllele : s.split(";")) {
        int delimiter = positionAllele.indexOf(':');
        long position = Long.parseLong(positionAllele.substring(0, delimiter));
        Integer idx = m_vcfPositionIndex.get(position);
        if (idx != null) {
          alleles[idx] = positionAllele.substring(delimiter + 1);
        }
      }
      return alleles;
    });
  }


  /**
   * Compares sample permutations with callable haplotypes using lazy BitSet intersections.
   * Match objects are accumulated by haplotype index and placed in a {@link TreeSet} only after their sequence sets
   * are complete because sequence content participates in {@code BaseMatch.compareTo()}.
   */
  protected SortedSet<HaplotypeMatch> comparePermutations() {
    initializeCandidateIndex();
    @Nullable HaplotypeMatch[] matches = new HaplotypeMatch[Objects.requireNonNull(m_haplotypeIndex).size()];
    for (SamplePermutation permutation : getPermutations()) {
      @Nullable String[] sequenceAlleles = permutation.getAllelesForMatching();
      BitSet candidates = new BitSet(m_haplotypeIndex.size());
      candidates.set(0, m_haplotypeIndex.size());
      for (int x = 0; x < sequenceAlleles.length; x += 1) {
        BitSet compatibleHaplotypes = getCompatibleHaplotypes(x, sequenceAlleles[x]);
        candidates.and(compatibleHaplotypes);
        if (candidates.isEmpty()) {
          break;
        }
      }
      for (int x = candidates.nextSetBit(0); x >= 0; x = candidates.nextSetBit(x + 1)) {
        if (matches[x] == null) {
          matches[x] = new HaplotypeMatch(getOutputHaplotype(x));
        }
        //noinspection DataFlowIssue
        matches[x].addSequence(permutation);
      }
    }
    return Arrays.stream(matches)
        .filter(Objects::nonNull)
        .collect(Collectors.toCollection(TreeSet::new));
  }


  /**
   * Initializes arrays aligned by haplotype index and empty per-position compatibility caches.
   * Expected alleles are resolved into permutation position order here so the hot loop can use array indexes only.
   */
  private void initializeCandidateIndex() {
    if (m_haplotypeIndex != null) {
      return;
    }
    m_haplotypeIndex = new ArrayList<>(getHaplotypes());
    m_haplotypeAlleles = new ArrayList<>(m_haplotypeIndex.size());
    @Nullable String[] referenceAlleles = null;
    if (m_defaultMissingAllelesToReference) {
      referenceAlleles = getReferenceHaplotype().getAlleles(m_permutationPositions);
    }
    for (NamedAllele haplotype : m_haplotypeIndex) {
      m_haplotypeAlleles.add(getAllelesForMatching(haplotype, m_permutationPositions, referenceAlleles));
    }
    m_outputHaplotypes = new ArrayList<>(Collections.nCopies(m_haplotypeIndex.size(), null));
    m_candidateIndex = new ArrayList<>(m_permutationPositions.length);
    for (int x = 0; x < m_permutationPositions.length; x += 1) {
      m_candidateIndex.add(new HashMap<>());
    }
  }


  private @Nullable String[] getAllelesForMatching(NamedAllele haplotype, VariantLocus[] positions,
      @Nullable String[] referenceAlleles) {
    @Nullable String[] alleles = haplotype.getAlleles(positions);
    if (!m_defaultMissingAllelesToReference || haplotype.isReference()) {
      return alleles;
    }
    assert referenceAlleles != null;
    return applyReferenceDefaults(positions, alleles, referenceAlleles);
  }


  private NamedAllele getOutputHaplotype(int index) {
    assert m_haplotypeIndex != null;
    NamedAllele haplotype = m_haplotypeIndex.get(index);
    if (!m_defaultMissingAllelesToReference || haplotype.isReference()) {
      return haplotype;
    }

    assert m_outputHaplotypes != null;
    NamedAllele outputHaplotype = m_outputHaplotypes.get(index);
    if (outputHaplotype == null) {
      outputHaplotype = materializeDefaultedHaplotype(haplotype, getReferenceHaplotype());
      m_outputHaplotypes.set(index, outputHaplotype);
    }
    return outputHaplotype;
  }


  private NamedAllele materializeDefaultedHaplotype(NamedAllele haplotype, NamedAllele referenceHaplotype) {
    @Nullable String[] alleles = applyReferenceDefaults(m_positions, haplotype.getAlleles(m_positions),
        referenceHaplotype.getAlleles(m_positions));
    @Nullable String[] cpicAlleles = applyReferenceCpicDefaults(haplotype, referenceHaplotype);
    NamedAllele outputHaplotype = new NamedAllele(haplotype.getId(), haplotype.getName(), alleles, cpicAlleles,
        haplotype.getMissingPositions(), haplotype.isReference());
    outputHaplotype.initialize(m_positions, haplotype.getScore());
    return outputHaplotype;
  }


  private NamedAllele getReferenceHaplotype() {
    return getHaplotypes().stream().filter(NamedAllele::isReference).findAny()
        .orElseThrow(() -> new IllegalStateException(m_gene + " does not have a reference"));
  }


  private @Nullable String[] applyReferenceDefaults(VariantLocus[] positions, @Nullable String[] alleles,
      @Nullable String[] referenceAlleles) {

    @Nullable String[] defaultedAlleles = new String[alleles.length];
    for (int x = 0; x < alleles.length; x += 1) {
      if (alleles[x] == null) {
        String refAllele = referenceAlleles[x];
        if (Iupac.isWobble(refAllele)) {
          defaultedAlleles[x] = positions[x].getRef();
        } else {
          defaultedAlleles[x] = refAllele;
        }
      } else {
        defaultedAlleles[x] = alleles[x];
      }
    }
    return defaultedAlleles;
  }


  private @Nullable String[] applyReferenceCpicDefaults(NamedAllele haplotype,
      NamedAllele referenceHaplotype) {

    @Nullable String[] cpicAlleles = new String[m_positions.length];
    for (int x = 0; x < m_positions.length; x += 1) {
      if (haplotype.getAllele(m_positions[x]) == null) {
        cpicAlleles[x] = referenceHaplotype.getCpicAlleles()[x];
      } else {
        cpicAlleles[x] = haplotype.getCpicAlleles()[x];
      }
    }
    return cpicAlleles;
  }


  /**
   * Gets haplotypes compatible with one observed allele, computing that cache entry on first use.
   * The returned BitSet is owned by the index; callers must clone it or use it only as the right-hand operand of a
   * BitSet operation.
   */
  private BitSet getCompatibleHaplotypes(int positionIndex, @Nullable String observedAllele) {
    assert m_candidateIndex != null;
    assert m_haplotypeIndex != null;
    assert m_haplotypeAlleles != null;
    Map<String, BitSet> candidatesByAllele = m_candidateIndex.get(positionIndex);
    return candidatesByAllele.computeIfAbsent(observedAllele, allele -> {
      BitSet compatibleHaplotypes = new BitSet(m_haplotypeIndex.size());
      for (int x = 0; x < m_haplotypeIndex.size(); x += 1) {
        if (m_haplotypeIndex.get(x).matchesAllele(m_haplotypeAlleles.get(x)[positionIndex], allele)) {
          compatibleHaplotypes.set(x);
        }
      }
      return compatibleHaplotypes;
    });
  }


  private void clearCandidateIndex() {
    m_haplotypeIndex = null;
    m_haplotypeAlleles = null;
    m_outputHaplotypes = null;
    m_candidateIndex = null;
  }


  @Override
  public String toString() {
    return m_gene + " match data for " + m_sampleId;
  }


  private static boolean areSampleAllelesHaploid(Collection<SampleAllele> sampleAlleles) {
    int s1 = 0;
    int s2 = 0;
    for (SampleAllele a : sampleAlleles) {
      if (a.getAllele1() == null) {
        s1 += 1;
      }
      if (a.getAllele2() == null) {
        s2 += 1;
      }
    }
    return s1 == sampleAlleles.size() || s2 == sampleAlleles.size();
  }
}
