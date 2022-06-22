package org.pharmgkb.pharmcat.util;

import java.util.Comparator;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


public class HaplotypeActivityComparator implements Comparator<Haplotype> {
  private static final Comparator<Haplotype> sf_comparator = new HaplotypeActivityComparator();

  public static Comparator<Haplotype> getComparator() {
    return sf_comparator;
  }

  @Override
  public int compare(Haplotype o1, Haplotype o2) {

    if (o1 == o2) {
      return 0;
    }
    if (o1 == null) {
      return -1;
    } else if (o2 == null) {
      return 1;
    }

    int rez = ObjectUtils.compare(o1.getGene(), o2.getGene());
    if (rez != 0) {
      return rez;
    }

    rez = ObjectUtils.compare(o1.getActivityValue(), o2.getActivityValue());
    if (rez != 0) {
      return rez;
    }

    return HaplotypeNameComparator.getComparator().compare(o1.getName(), o2.getName());
  }
}
