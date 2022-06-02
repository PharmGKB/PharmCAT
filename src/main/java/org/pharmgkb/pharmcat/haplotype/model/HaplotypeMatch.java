package org.pharmgkb.pharmcat.haplotype.model;

import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;


/**
 * This represents a haplotype and the sequences that matched it.
 *
 * @author Mark Woon
 */
public class HaplotypeMatch extends BaseMatch {


  public HaplotypeMatch(NamedAllele haplotype) {
    setName(haplotype.getName());
    setHaplotype(haplotype);
  }


  public boolean match(String seq) {
    if (getHaplotype().getPermutations().matcher(seq).matches()) {
      addSequence(seq);
      return true;
    }
    return false;
  }


  @Override
  public int compareTo(BaseMatch o) {
    if (this == o) {
      return 0;
    }
    if (o instanceof HaplotypeMatch hm) {
      int rez = HaplotypeNameComparator.getComparator().compare(getName(), hm.getName());
      if (rez != 0) {
        return rez;
      }
      rez = ObjectUtils.compare(getHaplotype(), hm.getHaplotype());
      if (rez != 0) {
        return rez;
      }
      rez = ObjectUtils.compare(getHaplotype().getScore(), hm.getHaplotype().getScore());
      if (rez != 0) {
        return rez;
      }
      return compareSequences(hm);
    }
    if (o instanceof CombinationMatch cm) {
      if (cm.getName().startsWith("g.")) {
        // push off-reference partial to bottom
        return -1;
      }
      int rez = ObjectUtils.compare(getHaplotype(), cm.getComponentHaplotypes().first());
      if (rez != 0) {
        return rez;
      }
      return -1;
    }
    throw new UnsupportedOperationException("Don't know how to compare HaplotypeMatch to " + o.getClass());
  }
}
