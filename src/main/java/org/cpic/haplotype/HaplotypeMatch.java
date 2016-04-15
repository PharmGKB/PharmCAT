package org.cpic.haplotype;

import java.util.HashSet;
import java.util.Set;
import javax.annotation.Nonnull;
import org.apache.commons.lang3.ObjectUtils;


/**
 * @author Mark Woon
 */
public class HaplotypeMatch implements Comparable<HaplotypeMatch> {
  private Haplotype m_haplotype;
  private Set<String> m_sequences = new HashSet<>();


  public HaplotypeMatch(@Nonnull Haplotype haplotype) {
    m_haplotype = haplotype;
  }

  public @Nonnull Haplotype getHaplotype() {
    return m_haplotype;
  }

  public boolean match(@Nonnull String seq) {
    if (m_haplotype.getPermutations().matcher(seq).matches()) {
      m_sequences.add(seq);
      return true;
    }
    return false;
  }

  public @Nonnull Set<String> getSequences() {
    return m_sequences;
  }


  @Override
  public int compareTo(@Nonnull HaplotypeMatch o) {

    int rez = ObjectUtils.compare(m_haplotype, o.getHaplotype());
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_sequences.size(), o.getSequences().size());
  }

  @Override
  public String toString() {
    return m_haplotype.toString();
  }
}
