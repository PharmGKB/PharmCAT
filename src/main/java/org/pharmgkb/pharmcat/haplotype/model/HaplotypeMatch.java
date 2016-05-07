package org.pharmgkb.pharmcat.haplotype.model;

import java.util.HashSet;
import java.util.Set;
import javax.annotation.Nonnull;
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
  @SerializedName("name")
  @Expose
  private String m_name;
  private NamedAllele m_haplotype;
  @Expose
  private Set<String> m_sequences = new HashSet<>();


  public HaplotypeMatch(@Nonnull NamedAllele haplotype) {
    m_haplotype = haplotype;
    m_name = m_haplotype.getName();
  }


  public String getName() {
    return m_name;
  }

  public @Nonnull NamedAllele getHaplotype() {
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
