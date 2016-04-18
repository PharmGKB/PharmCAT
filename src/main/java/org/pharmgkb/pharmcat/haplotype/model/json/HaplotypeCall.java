package org.pharmgkb.pharmcat.haplotype.model.json;

import java.util.Set;


/**
 * A haplotype call result.
 *
 * @author Mark Woon
 */
public class HaplotypeCall {
  private String m_name;
  private Set<String> m_alleles;


  public HaplotypeCall(String name, Set<String> alleles) {
    m_name = name;
    m_alleles = alleles;
  }

  public String getName() {
    return m_name;
  }


  /**
   * Gets the allele sequences that matched this haplotype.
   * This is mainly useful when dealing with unphased data.
   */
  public Set<String> getAlleles() {
    return m_alleles;
  }
}
