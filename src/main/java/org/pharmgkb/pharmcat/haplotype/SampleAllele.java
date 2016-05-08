package org.pharmgkb.pharmcat.haplotype;

import javax.annotation.Nonnull;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.common.comparator.ChromosomeNameComparator;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.definition.model.VariantType;


/**
 * This class represents a sample's allele at a specific location.
 * <p>
 * It's comparator only takes it's chromosomal position into account.
 *
 * @author Mark Woon
 */
public class SampleAllele implements Comparable<SampleAllele> {
  private String m_chromosome;
  private int m_position;
  private String m_allele1;
  private String m_allele2;
  private boolean m_isPhased;

  public SampleAllele(@Nonnull String chromosome, long position, @Nonnull String a1, String a2, boolean isPhased) {
    m_chromosome = chromosome;
    m_position = (int)position;
    m_allele1 = a1;
    m_allele2 = a2;
    m_isPhased = isPhased;
  }

  public String getChromosome() {
    return m_chromosome;
  }

  public int getPosition() {
    return m_position;
  }

  public String getChrPosition() {
    return m_chromosome + ":" + m_position;
  }

  public String getAllele1() {
    return m_allele1;
  }

  public String getAllele2() {
    return m_allele2;
  }

  public boolean isPhased() {
    return m_isPhased;
  }

  @Override
  public String toString() {
    return m_allele1 + "/" + m_allele2 + " @ " + m_chromosome + ":" + m_position;
  }

  @Override
  public int compareTo(@Nonnull SampleAllele o) {

    int rez = ChromosomeNameComparator.getComparator().compare(m_chromosome, o.getChromosome());
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_position, o.getPosition());
  }


  /**
   * Interprets the alleles in this {@link SampleAllele} in terms of the given {@link VariantLocus}.
   * This will return a <strong>new</strong> {@link SampleAllele} if the {@link VariantLocus} is not a SNP, with it's
   * alleles modified to use the format used by the allele definitions.
   */
  public SampleAllele forVariant(VariantLocus variant) {

    if (variant.getType() == VariantType.SNP) {
      return this;
    }
    String a1 = m_allele1;
    String a2 = m_allele2;

    if (variant.getType() == VariantType.INS) {
      // VCF:         TC  -> TCA
      // definition:  del -> insA
      if (a1.length() > a2.length()) {
        a1 = "ins" + a1.substring(a2.length());
        a2 = "del";
      } else {
        a2 = "ins" + a2.substring(a1.length());
        a1 = "del";
      }

    } else if (variant.getType() == VariantType.DEL) {
      // VCF:         TC -> T
      // definition:  C  -> delC
      if (a1.length() > a2.length()) {
        a1 = a1.substring(1);
        a2 = "del" + a1;
      } else {
        a2 = a2.substring(1);
        a1 = "del" + a2;
      }
    }
    return new SampleAllele(m_chromosome, m_position, a1, a2, m_isPhased);
  }
}
