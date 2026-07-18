package org.pharmgkb.pharmcat.haplotype;

import java.util.Arrays;
import org.jspecify.annotations.Nullable;


/**
 * One possible single-strand interpretation of the sample genotype.
 *
 * <p>Positions and alleles are parallel arrays in numeric position order. Matching uses the arrays directly; the
 * encoded sequence is created lazily only for stable result ordering and compatibility with legacy match paths.</p>
 */
public final class SamplePermutation {
  private final long[] m_positions;
  private final @Nullable String[] m_alleles;
  private @Nullable String m_sequence;


  SamplePermutation(@Nullable String[] alleles, long[] positions) {
    m_alleles = Arrays.copyOf(alleles, alleles.length);
    m_positions = positions;
  }


  /**
   * Gets a copy of alleles in position order.
   * The internal array must remain immutable because it is used for equality, hashing, and lazy sequence generation.
   */
  public @Nullable String[] getAlleles() {
    return Arrays.copyOf(m_alleles, m_alleles.length);
  }

  /**
   * Package-private raw array access for matcher hot paths.
   * Do not expose this outside the haplotype package; callers must not mutate the returned array.
   */
  @Nullable String[] getAllelesForMatching() {
    return m_alleles;
  }


  public String getSequence() {
    if (m_sequence == null) {
      StringBuilder builder = new StringBuilder();
      for (int x = 0; x < m_positions.length; x += 1) {
        if (!builder.isEmpty()) {
          builder.append(";");
        }
        builder.append(m_positions[x])
            .append(":")
            .append(m_alleles[x]);
      }
      m_sequence = builder.toString();
    }
    return m_sequence;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (!(o instanceof SamplePermutation that)) {
      return false;
    }
    return Arrays.equals(m_positions, that.m_positions) && Arrays.equals(m_alleles, that.m_alleles);
  }


  @Override
  public int hashCode() {
    return 31 * Arrays.hashCode(m_positions) + Arrays.hashCode(m_alleles);
  }
}
