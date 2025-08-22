package org.pharmgkb.pharmcat.definition.model;

import java.util.SortedSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * This class represents an unphased diplotype priority override.
 *
 * @author Mark Woon
 */
public class UnphasedDiplotypePriority implements Comparable<UnphasedDiplotypePriority> {
  @Expose
  @SerializedName("id")
  private final String m_id;
  @Expose
  @SerializedName("list")
  private final SortedSet<String> m_list;
  @Expose
  @SerializedName("pick")
  private final String m_pick;


  UnphasedDiplotypePriority(SortedSet<String> list, String pick) {
    if (!list.contains(pick)) {
      throw new IllegalArgumentException("'" + pick + "' is not contained within list of diplotypes (" +
          String.join(",", list) + ")");
    }
    m_id = DefinitionExemption.generateUnphasedPriorityKey(list);
    this.m_list = list;
    this.m_pick = pick;
  }

  public String getId() {
    return m_id;
  }

  public String getPick() {
    return m_pick;
  }

  public SortedSet<String> getList() {
    return m_list;
  }

  @Override
  public int compareTo(UnphasedDiplotypePriority o) {
    if (this == o) {
      return 0;
    }
    String[] haps1 = m_pick.split("/");
    String[] haps2 = o.m_pick.split("/");

    int rez = HaplotypeNameComparator.getComparator().compare(haps1[0], haps2[0]);
    if (rez != 0) {
      return rez;
    }
    rez = HaplotypeNameComparator.getComparator().compare(haps1[1], haps2[1]);
    if (rez != 0) {
      return rez;
    }

    return m_id.compareTo(o.m_id);
  }
}
