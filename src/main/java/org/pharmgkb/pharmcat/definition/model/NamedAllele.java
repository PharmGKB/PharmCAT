package org.pharmgkb.pharmcat.definition.model;

import java.util.Map;

/**
 * A named allele (aka Haplotype, Star Allele, etc.).
 *
 * @author Ryan Whaley
 */
public class NamedAllele {

  private String m_name;
  private String m_id;
  private String m_function;
  private String[] m_alleles;
  private Map<String,String> m_popFreqMap;

  public String toString() {
    return m_name + " [" + m_id + "]";
  }

  /**
   * The function of this allele (e.g. normal, reduced function)
   */
  public String getFunction() {
    return m_function;
  }

  public void setFunction(String function) {
    m_function = function;
  }

  /**
   * The name of this named allele (e.g. *1, Foo123Bar)
   */
  public String getName() {
    return m_name;
  }

  public void setName(String name) {
    m_name = name;
  }

  /**
   * The CPIC identifier for this named allele (e.g. CA10000.1)
   */
  public String getId() {
    return m_id;
  }

  public void setId(String id) {
    m_id = id;
  }

  /**
   * The array of alleles that define this allele.
   *
   * <em>Note:</em> use this in conjunction with <code>AlleleTranslation.variants</code> to get the name of the variant
   */
  public String[] getAlleles() {
    return m_alleles;
  }

  public void setAlleles(String[] alleles) {
    m_alleles = alleles;
  }

  /**
   * A mapping of population name to allele frequency
   */
  public Map<String, String> getPopFreqMap() {
    return m_popFreqMap;
  }

  public void setPopFreqMap(Map<String, String> popFreqMap) {
    m_popFreqMap = popFreqMap;
  }
}
