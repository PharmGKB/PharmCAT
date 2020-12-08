package org.pharmgkb.pharmcat.definition.model;

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
  private String m_generesult;
  @SerializedName("diplotype")
  @Expose
  private String m_diplotype;
  @SerializedName("description")
  @Expose
  private String m_description;
  @SerializedName("lookupkey")
  @Expose
  private String m_lookupKey;

  public String getGeneresult() {
    return m_generesult;
  }

  public void setGeneresult(String generesult) {
    m_generesult = generesult;
  }

  /**
   * The bare diplotype name
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
}
