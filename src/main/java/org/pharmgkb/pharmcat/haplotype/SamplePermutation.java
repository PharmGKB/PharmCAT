package org.pharmgkb.pharmcat.haplotype;

import java.util.Arrays;
import org.jspecify.annotations.Nullable;


/**
 * A sample haplotype permutation represented by its alleles at each position.
 */
final class SamplePermutation {
  private final long[] m_positions;
  private final @Nullable String[] m_alleles;
  private @Nullable String m_sequence;


  SamplePermutation(@Nullable String[] alleles, long[] positions) {
    m_alleles = Arrays.copyOf(alleles, alleles.length);
    m_positions = positions;
  }


  @Nullable String[] getAlleles() {
    return m_alleles;
  }


  String getSequence() {
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
