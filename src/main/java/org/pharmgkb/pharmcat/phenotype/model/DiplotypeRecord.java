package org.pharmgkb.pharmcat.phenotype.model;

import java.util.Map;
import com.google.common.base.Objects;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * Model that represents a diplotype data in the CPIC database
 *
 * @author Ryan Whaley
 */
public class DiplotypeRecord {

  @SerializedName("generesult")
  @Expose
  private String m_geneResult;
  @SerializedName("diplotype")
  @Expose
  private String m_diplotype;
  @SerializedName("description")
  @Expose
  private String m_description;
  @SerializedName("lookupkey")
  @Expose
  private String m_lookupKey;
  @SerializedName("diplotypekey")
  @Expose
  private Map<String,Integer> m_diplotypeKey;


  public String getGeneResult() {
    return m_geneResult;
  }

  public void setGeneResult(String geneResult) {
    m_geneResult = geneResult;
  }

  /**
   * The bare diplotype name.
   *
   * @return diplotype in the form "*1/*3"
   */
  public String getDiplotype() {
    return m_diplotype;
  }

  public void setDiplotype(String diplotype) {
    m_diplotype = diplotype;
  }

  public String toString() {
    return m_diplotype;
  }

  public String getDescription() {
    return m_description;
  }

  public void setDescription(String description) {
    m_description = description;
  }

  public String getLookupKey() {
    return m_lookupKey;
  }

  public void setLookupKey(String lookupKey) {
    m_lookupKey = lookupKey;
  }

  public Map<String, Integer> getDiplotypeKey() {
    return m_diplotypeKey;
  }

  public void setDiplotypeKey(Map<String, Integer> diplotypeKey) {
    m_diplotypeKey = diplotypeKey;
  }


  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj == this) {
      return true;
    }
    if (obj.getClass() != getClass()) {
      return false;
    }
    DiplotypeRecord other = (DiplotypeRecord)obj;
    return Objects.equal(m_geneResult, other.getGeneResult()) &&
        Objects.equal(m_diplotype, other.getDiplotype()) &&
        Objects.equal(m_description, other.getDescription()) &&
        Objects.equal(m_lookupKey, other.getLookupKey()) &&
        Objects.equal(m_diplotypeKey, other.getDiplotypeKey());
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(m_geneResult, m_diplotype, m_description, m_lookupKey, m_diplotypeKey);
  }
}
