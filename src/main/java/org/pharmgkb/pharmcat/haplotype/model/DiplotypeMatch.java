package org.pharmgkb.pharmcat.haplotype.model;

import java.util.HashSet;
import java.util.Set;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.ObjectUtils;


/**
 * This represents a diplotype and the sequences that matched it.
 *
 * @author Mark Woon
 */
public class DiplotypeMatch implements Comparable<DiplotypeMatch> {
  private String m_name;
  private HaplotypeMatch m_haplotype1;
  private HaplotypeMatch m_haplotype2;
  private int m_score;
  private Set<String[]> m_sequences = new HashSet<>();



  public DiplotypeMatch(@Nonnull HaplotypeMatch hm1, @Nonnull HaplotypeMatch hm2) {
    m_haplotype1 = hm1;
    m_haplotype2 = hm2;
    m_name = m_haplotype1.getName() + "/" + m_haplotype2.getName();
    m_score = m_haplotype1.getHaplotype().getVariants().size() + m_haplotype2.getHaplotype().getVariants().size();
  }

  public String getName() {
    return m_name;
  }

  public int getScore() {
    return m_score;
  }

  public HaplotypeMatch getHaplotype1() {
    return m_haplotype1;
  }

  public HaplotypeMatch getHaplotype2() {
    return m_haplotype2;
  }

  public Set<String[]> getSequences() {
    return m_sequences;
  }

  public void addSequencePair(@Nonnull  String[] pair) {
    Preconditions.checkNotNull(pair);
    Preconditions.checkArgument(pair.length == 2, "Sequence pair must have 2 sequences");
    m_sequences.add(pair);
  }

  @Override
  public String toString() {
    return m_name;
  }

  @Override
  public int compareTo(@Nonnull DiplotypeMatch o) {

    int rez = ObjectUtils.compare(m_score, o.getScore());
    if (rez != 0) {
      return rez;
    }
    rez = ObjectUtils.compare(m_haplotype1, o.getHaplotype1());
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_haplotype2, o.getHaplotype2());
  }
}
