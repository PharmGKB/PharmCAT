package org.pharmgkb.pharmcat.definition.model;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.haplotype.Iupac;
import org.pharmgkb.pharmcat.reporter.TextConstants;


/**
 * Immutable per-gene BitSet index used by {@code MatchData.comparePermutations()} to look up haplotypes compatible with
 * an observed allele in constant time.
 *
 * <p>The index depends only on the callable {@link NamedAllele} set, the ordered permutation positions, and the
 * "default missing alleles to reference" flag; when the sample has no missing positions, all three are determined by
 * (gene, findCombinations, defaultMissingAllelesToReference), which lets {@link DefinitionFile} cache and share
 * instances across samples.</p>
 *
 * <p>Fields correspond to the three-way partition originally documented on {@code MatchData}:</p>
 * <ul>
 *   <li>{@link #wildcardsAt} - haplotypes whose raw allele is {@code null} and matches any observed value.</li>
 *   <li>{@link #defaultedToRefAt} - non-reference haplotypes whose null slot defaults to reference; matched only when
 *       the observed allele expands to the reference-defaulted value.</li>
 *   <li>{@link #specificsAt} - haplotypes with a non-null raw allele, wobble-expanded at build time.</li>
 * </ul>
 */
public final class HaplotypeCandidateIndex {
  private final List<NamedAllele> m_haplotypeIndex;
  private final VariantLocus[] m_permutationPositions;
  private final BitSet[] m_wildcardsAt;
  private final BitSet[] m_defaultedToRefAt;
  private final Set<String>[] m_refExpandedBasesAt;
  private final List<Map<String, BitSet>> m_specificsAt;
  private final boolean m_defaultMissingAllelesToReference;


  public HaplotypeCandidateIndex(SortedSet<NamedAllele> haplotypes, VariantLocus[] permutationPositions,
      boolean defaultMissingAllelesToReference) {
    m_haplotypeIndex = List.copyOf(haplotypes);
    m_permutationPositions = permutationPositions;
    m_defaultMissingAllelesToReference = defaultMissingAllelesToReference;

    int numHaps = m_haplotypeIndex.size();
    int numPos = permutationPositions.length;

    m_wildcardsAt = new BitSet[numPos];
    m_defaultedToRefAt = new BitSet[numPos];
    m_refExpandedBasesAt = new Set[numPos];
    m_specificsAt = new ArrayList<>(numPos);
    for (int p = 0; p < numPos; p += 1) {
      m_wildcardsAt[p] = new BitSet(numHaps);
      m_defaultedToRefAt[p] = new BitSet(numHaps);
      m_refExpandedBasesAt[p] = Collections.emptySet();
      m_specificsAt.add(new HashMap<>());
    }

    @Nullable String[] refExpected = null;
    if (defaultMissingAllelesToReference) {
      NamedAllele referenceHaplotype = m_haplotypeIndex.stream()
          .filter(NamedAllele::isReference)
          .findAny()
          .orElseThrow(() -> new IllegalStateException("No reference haplotype in candidate set"));
      refExpected = referenceHaplotype.getAlleles(permutationPositions);
      for (int p = 0; p < numPos; p += 1) {
        String refA = refExpected[p];
        if (refA == null) {
          continue;
        }
        String refValue = Iupac.isWobble(refA) ? permutationPositions[p].getRef() : refA;
        m_refExpandedBasesAt[p] = new HashSet<>(expandAllele(refValue));
      }
    }

    for (int h = 0; h < numHaps; h += 1) {
      NamedAllele hap = m_haplotypeIndex.get(h);
      @Nullable String[] rawAlleles = hap.getAlleles(permutationPositions);
      boolean isRef = hap.isReference();
      for (int p = 0; p < numPos; p += 1) {
        String raw = rawAlleles[p];
        if (raw != null) {
          for (String base : expandAllele(raw)) {
            m_specificsAt.get(p).computeIfAbsent(base, k -> new BitSet(numHaps)).set(h);
          }
        } else if (defaultMissingAllelesToReference && !isRef && refExpected != null && refExpected[p] != null) {
          m_defaultedToRefAt[p].set(h);
        } else {
          m_wildcardsAt[p].set(h);
        }
      }
    }
  }


  public int numHaplotypes() {
    return m_haplotypeIndex.size();
  }

  public int numPositions() {
    return m_permutationPositions.length;
  }

  public NamedAllele haplotypeAt(int index) {
    return m_haplotypeIndex.get(index);
  }

  public List<NamedAllele> haplotypes() {
    return m_haplotypeIndex;
  }

  public BitSet wildcardsAt(int p) {
    return m_wildcardsAt[p];
  }

  public BitSet defaultedToRefAt(int p) {
    return m_defaultedToRefAt[p];
  }

  public Set<String> refExpandedBasesAt(int p) {
    return m_refExpandedBasesAt[p];
  }

  public Map<String, BitSet> specificsAt(int p) {
    return m_specificsAt.get(p);
  }

  public boolean isDefaultMissingAllelesToReference() {
    return m_defaultMissingAllelesToReference;
  }


  /**
   * Expands a raw haplotype allele into the set of observed-allele strings it matches, mirroring
   * {@link NamedAllele#matchesAllele} semantics exactly.
   */
  private static List<String> expandAllele(String raw) {
    if (raw.contains(TextConstants.REPEAT_WOBBLE_DELIMITER)) {
      return Arrays.asList(raw.split(TextConstants.REPEAT_WOBBLE_DELIMITER));
    }
    if (raw.length() == 1) {
      Iupac iupac = Iupac.lookup(raw);
      if (iupac.isAmbiguity()) {
        return iupac.getBases();
      }
      return List.of(iupac.getRegex());
    }
    return List.of(raw);
  }
}
