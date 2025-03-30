package org.pharmgkb.pharmcat.haplotype.model;

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
}
