package org.pharmgkb.pharmcat.haplotype.model;

import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;


/**
 * This represents a haplotype and the sequences that matched it.
 *
 * @author Mark Woon
 */
public class HaplotypeMatch extends BaseMatch {
  private final @Nullable String @Nullable[] m_expectedAlleles;


  public HaplotypeMatch(NamedAllele haplotype) {
    setName(haplotype.getName());
    setHaplotype(haplotype);
    m_expectedAlleles = null;
  }


  public HaplotypeMatch(NamedAllele haplotype, VariantLocus[] positions) {
    setName(haplotype.getName());
    setHaplotype(haplotype);
    m_expectedAlleles = haplotype.getAlleles(positions);
  }


  /**
   * Matches an encoded sequence using the legacy regex representation.
   */
  @Deprecated
  public boolean match(String seq) {
    if (getHaplotype().getPermutations().matcher(seq).matches()) {
      addSequence(seq);
      return true;
    }
    return false;
  }


  public boolean match(String seq, VariantLocus[] positions, String[] sequenceAlleles) {
    @Nullable String[] expectedAlleles = m_expectedAlleles;
    if (expectedAlleles == null) {
      expectedAlleles = getHaplotype().getAlleles(positions);
    }
    if (getHaplotype().matches(expectedAlleles, sequenceAlleles)) {
      addSequence(seq);
      return true;
    }
    return false;
  }
}
