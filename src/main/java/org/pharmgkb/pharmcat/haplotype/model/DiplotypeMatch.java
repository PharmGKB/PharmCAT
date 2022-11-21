package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Objects;
import java.util.Set;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
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
  private final BaseMatch m_haplotype1;
  @Expose
  @SerializedName("haplotype2")
  private BaseMatch m_haplotype2;
  @Expose
  @SerializedName("score")
  private int m_score;
  private final Set<String[]> m_sequences = new HashSet<>();
  private final MatchData m_dataset;



  public DiplotypeMatch(BaseMatch hm1, @Nullable BaseMatch hm2, MatchData dataset) {
    if (hm2 != null) {
      BaseMatch[] sortedHaps = new BaseMatch[]{ hm1, hm2 };
      Arrays.sort(sortedHaps);
      m_haplotype1 = sortedHaps[0];
      m_haplotype2 = sortedHaps[1];
      m_name = m_haplotype1.getName() + "/" + m_haplotype2.getName();
      m_score = m_haplotype1.getHaplotype().getScore() + m_haplotype2.getHaplotype().getScore();
    } else {
      // haploid
      m_haplotype1 = hm1;
      m_name = m_haplotype1.getName();
      m_score = m_haplotype1.getHaplotype().getScore();
    }
    m_dataset = dataset;
  }

  public String getName() {
    return m_name;
  }

  public int getScore() {
    return m_score;
  }

  public void setScore(int score) {
    m_score = score;
  }

  public BaseMatch getHaplotype1() {
    return m_haplotype1;
  }

  public @Nullable BaseMatch getHaplotype2() {
    return m_haplotype2;
  }

  public Set<String[]> getSequences() {
    return m_sequences;
  }

  public void addSequencePair(String[] pair) {
    Preconditions.checkNotNull(pair);
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
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj == this) {
      return true;
    }
    if (obj.getClass() != getClass()) {
      return false;
    }
    DiplotypeMatch dm = (DiplotypeMatch)obj;
    return Objects.equals(m_score, dm.getScore()) &&
        Objects.equals(m_haplotype1, dm.getHaplotype1()) &&
        Objects.equals(m_haplotype2, dm.getHaplotype2());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_score, m_haplotype1, m_haplotype2);
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
