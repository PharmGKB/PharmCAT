package org.pharmgkb.pharmcat.haplotype.model;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.haplotype.MatchData;


/**
 * This represents a diplotype and the sequences that matched it.
 *
 * @author Mark Woon
 */
public class DiplotypeMatch implements Comparable<DiplotypeMatch> {
  @Expose
  @SerializedName("name")
  private final String m_name;
  @Expose
  @SerializedName("haplotype1")
  private final HaplotypeMatch m_haplotype1;
  @Expose
  @SerializedName("haplotype2")
  private final HaplotypeMatch m_haplotype2;
  @Expose
  @SerializedName("score")
  private final int m_score;
  private final Set<String[]> m_sequences = new HashSet<>();
  private final MatchData m_dataset;



  public DiplotypeMatch(HaplotypeMatch hm1, HaplotypeMatch hm2, MatchData dataset) {
    m_haplotype1 = hm1;
    m_haplotype2 = hm2;
    List<String> hapNames = new ArrayList<>();
    hapNames.add(m_haplotype1.getName());
    hapNames.add(m_haplotype2.getName());
    hapNames.sort(HaplotypeNameComparator.getComparator());
    m_name = String.join("/", hapNames);
    m_score = m_haplotype1.getHaplotype().getScore() + m_haplotype2.getHaplotype().getScore();
    m_dataset = dataset;
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

  public void addSequencePair(String[] pair) {
    Preconditions.checkNotNull(pair);
    Preconditions.checkArgument(pair.length == 2, "Sequence pair must have 2 sequences");
    m_sequences.add(pair);
  }

  public MatchData getDataset() {
    return m_dataset;
  }


  @Override
  public String toString() {
    return m_name;
  }

  @Override
  public int compareTo(DiplotypeMatch o) {

    int rez = ObjectUtils.compare(o.getScore(), m_score);
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
