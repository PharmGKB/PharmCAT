package org.pharmgkb.pharmcat.haplotype;

import java.util.List;
import org.apache.commons.lang3.ObjectUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.comparator.ChromosomeNameComparator;


/**
 * This class represents a sample's allele at a specific location.
 * <p>
 * It's comparator only takes it's chromosomal position into account.
 *
 * @author Mark Woon
 */
public class SampleAllele implements Comparable<SampleAllele> {
  private final String m_chromosome;
  private final int m_position;
  private final String m_allele1;
  private String m_allele2;
  private final boolean m_isPhased;
  private final List<String> m_vcfAlleles;


  public SampleAllele(String chromosome, long position, String a1, @Nullable String a2,
      boolean isPhased, List<String> vcfAlleles) {
    m_chromosome = chromosome;
    m_position = (int)position;
    m_allele1 = a1.toUpperCase();
    if (a2 != null) {
      m_allele2 = a2.toUpperCase();
    }
    m_isPhased = isPhased;
    m_vcfAlleles = vcfAlleles;
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

  public List<String> getVcfAlleles() {
    return m_vcfAlleles;
  }

  @Override
  public String toString() {
    return m_allele1 + "/" + m_allele2 + " @ " + m_chromosome + ":" + m_position;
  }

  @Override
  public int compareTo(SampleAllele o) {

    int rez = ChromosomeNameComparator.getComparator().compare(m_chromosome, o.getChromosome());
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_position, o.getPosition());
  }
}
