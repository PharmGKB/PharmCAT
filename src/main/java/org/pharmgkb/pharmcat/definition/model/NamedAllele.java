package org.pharmgkb.pharmcat.definition.model;

import java.util.Arrays;
import java.util.Map;
import com.google.common.base.Objects;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.ObjectUtils;
import org.pharmgkb.common.comparator.HaplotypeNameComparator;


/**
 * A named allele (aka Haplotype, Star Allele, etc.).
 *
 * @author Ryan Whaley
 */
public class NamedAllele implements Comparable<NamedAllele> {
  @SerializedName("name")
  private String m_name;
  @SerializedName("id")
  private String m_id;
  @SerializedName("function")
  private String m_function;
  @SerializedName("alleles")
  private String[] m_alleles;
  @SerializedName("populationFrequency")
  private Map<String,String> m_popFreqMap;


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
   * The function of this allele (e.g. normal, reduced function).
   */
  public String getFunction() {
    return m_function;
  }

  public void setFunction(String function) {
    m_function = function;
  }


  /**
   * The array of alleles that define this allele.
   *
   * <em>Note:</em> use this in conjunction with {@link DefinitionFile#getVariants()} to get the name of the variant
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


  @Override
  public String toString() {
    return m_name + " [" + m_id + "]";
  }

  @Override
  public int compareTo(NamedAllele o) {
    int rez = HaplotypeNameComparator.getComparator().compare(m_name, o.m_name);
    if (rez != 0) {
      return rez;
    }
    return ObjectUtils.compare(m_id, o.m_id);
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (!(o instanceof NamedAllele)) {
      return false;
    }
    NamedAllele that = (NamedAllele)o;
    return Objects.equal(m_name, that.getName()) &&
        Objects.equal(m_id, that.getId()) &&
        Objects.equal(m_function, that.getFunction()) &&
        Arrays.equals(m_alleles, that.getAlleles()) &&
        Objects.equal(m_popFreqMap, that.getPopFreqMap());
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(m_name, m_id, m_function, m_alleles, m_popFreqMap);
  }
}
