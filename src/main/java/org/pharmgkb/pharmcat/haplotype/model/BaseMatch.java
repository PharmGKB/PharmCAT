package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.collect.Lists;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * This is the base class for a {@link NamedAllele} match on a single strand.
 *
 * @author Mark Woon
 */
public class BaseMatch implements Comparable<BaseMatch> {
  @Expose
  @SerializedName("name")
  private String m_name;
  @Expose
  @SerializedName("haplotype")
  private NamedAllele m_haplotype;
  @Expose
  @SerializedName("sequences")
  private final SortedSet<String> m_sequences = new TreeSet<>();


  public String getName() {
    return m_name;
  }

  void setName(String name) {
    m_name = name;
  }


  public NamedAllele getHaplotype() {
    return m_haplotype;
  }

  void setHaplotype(NamedAllele haplotype) {
    m_haplotype = haplotype;
  }


  /**
   * Gets sequences of sample alleles that match this haplotype.
   * May have multiple matches due to unphased data.
   */
  public SortedSet<String> getSequences() {
    return m_sequences;
  }

  public void addSequence(String seq) {
    m_sequences.add(seq);
  }


  @Override
  public int compareTo(BaseMatch o) {
    if (this == o) {
      return 0;
    }
    String name1 = getName();
    String name2 = o.getName();

    int rez = HaplotypeNameComparator.getComparator().compare(name1, name2);
    if (rez != 0) {
      return rez;
    }

    if (this instanceof HaplotypeMatch && o instanceof HaplotypeMatch hm2) {
      rez = ObjectUtils.compare(getHaplotype().getScore(), hm2.getHaplotype().getScore());
      if (rez != 0) {
        return rez;
      }
    }

    return compareSequences(o);
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (!(obj instanceof BaseMatch that)) {
      return false;
    }
    return Objects.equals(getName(), that.getName()) &&
        compareSequences(that) == 0;
  }

  @Override
  public int hashCode() {
    return Objects.hash(getName(), Arrays.hashCode(getSequences().toArray()));
  }

  int compareSequences(BaseMatch o) {

    int rez = Integer.compare(getSequences().size(), o.getSequences().size());
    if (rez != 0) {
      return rez;
    }
    Iterator<String> it1 = getSequences().iterator();
    Iterator<String> it2 = o.getSequences().iterator();
    while (it1.hasNext()) {
      String n1 = it1.next();
      String n2 = it2.next();
      rez = n1.compareTo(n2);
      if (rez != 0) {
        return rez;
      }
    }
    return 0;
  }


  @Override
  public String toString() {
    return getName();
  }


  /**
   * Gets a list of all haplotype names involved in this match.
   */
  public List<String> getHaplotypeNames() {
    if (this instanceof CombinationMatch cm) {
      return cm.getComponentHaplotypes().stream()
          .map(NamedAllele::getName)
          .collect(Collectors.toList());
    } else {
      return Lists.newArrayList(getHaplotype().getName());
    }
  }
}
