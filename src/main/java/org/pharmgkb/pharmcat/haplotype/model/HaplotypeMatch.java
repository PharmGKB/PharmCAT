package org.pharmgkb.pharmcat.haplotype.model;

import java.util.SortedSet;
import java.util.TreeSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;


/**
 * This represents a haplotype and the sequences that matched it.
 *
 * @author Mark Woon
 */
public class HaplotypeMatch implements Comparable<HaplotypeMatch> {
  @Expose
  @SerializedName("name")
  private String m_name;
  private NamedAllele m_haplotype;
  @Expose
  @SerializedName("sequences")
  private SortedSet<String> m_sequences = new TreeSet<>();


  public HaplotypeMatch(NamedAllele haplotype) {
    m_haplotype = haplotype;
    m_name = m_haplotype.getName();
  }


  public String getName() {
    return m_name;
  }

  public NamedAllele getHaplotype() {
    return m_haplotype;
  }


  public boolean match(String seq) {
    if (m_haplotype.getPermutations().matcher(seq).matches()) {
      m_sequences.add(seq);
      return true;
    }
    return false;
  }

  public SortedSet<String> getSequences() {
    return m_sequences;
  }


  @Override
  public int compareTo(HaplotypeMatch o) {

    int rez = ObjectUtils.compare(m_name, o.getName());
    if (rez != 0) {
      return rez;
    }
    rez = ObjectUtils.compare(m_haplotype, o.getHaplotype());
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_sequences.size(), o.getSequences().size());
  }

  @Override
  public String toString() {
    return m_name;
  }
}
