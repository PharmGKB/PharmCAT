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


  public boolean matches(String[] sequenceAlleles) {
    if (m_expectedAlleles == null) {
      throw new IllegalStateException("No expected alleles were initialized for this match");
    }
    return getHaplotype().matches(m_expectedAlleles, sequenceAlleles);
  }
}
