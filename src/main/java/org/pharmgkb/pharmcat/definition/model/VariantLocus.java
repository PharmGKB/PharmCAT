package org.pharmgkb.pharmcat.definition.model;

import java.util.Objects;
import com.google.gson.annotations.SerializedName;


/**
 * All the information used to describe a particular location of a variant.
 *
 * @author Ryan Whaley
 */
public class VariantLocus {
  @SerializedName("rsid")
  private String m_rsid;
  @SerializedName("chrPosition")
  private String m_chrPosition;
  @SerializedName("genePosition")
  private String m_genePosition;
  @SerializedName("proteinNote")
  private String m_proteinNote;
  @SerializedName("resourceNote")
  private String m_resourceNote;
  @SerializedName("isInDel")
  private boolean m_isInDel;


  /**
   * A common name for this variant, usually specified by some specialized resource
   */
  public String getResourceNote() {
    return m_resourceNote;
  }

  public void setResourceNote(String resourceNote) {
    m_resourceNote = resourceNote;
  }

  /**
   * The name use for this location on the chromosomal sequence, should be relative to plus strand
   */
  public String getChrPosition() {
    return m_chrPosition;
  }

  public void setChrPosition(String chrPosition) {
    m_chrPosition = chrPosition;
  }

  /**
   * The name use for this location on the gene sequence, relative to the strand the gene is on
   */
  public String getGenePosition() {
    return m_genePosition;
  }

  public void setGenePosition(String genePosition) {
    m_genePosition = genePosition;
  }

  /**
   * The name use for this location on the protein sequence
   */
  public String getProteinNote() {
    return m_proteinNote;
  }

  public void setProteinNote(String proteinNote) {
    m_proteinNote = proteinNote;
  }

  /**
   * The identifier use for this location from dbSNP
   */
  public String getRsid() {
    return m_rsid;
  }

  public void setRsid(String rsid) {
    m_rsid = rsid;
  }


  public boolean isInDel() {
    return m_isInDel;
  }

  public void setInDel(boolean inDel) {
    m_isInDel = inDel;
  }


  @Override
  public boolean equals(Object o) {
    if (this == o) {
      return true;
    }
    if (!(o instanceof VariantLocus)) {
      return false;
    }
    VariantLocus that = (VariantLocus)o;
    return m_isInDel == that.isInDel() &&
        Objects.equals(m_rsid, that.getRsid()) &&
        Objects.equals(m_chrPosition, that.getChrPosition()) &&
        Objects.equals(m_genePosition, that.getGenePosition()) &&
        Objects.equals(m_proteinNote, that.getProteinNote()) &&
        Objects.equals(m_resourceNote, that.getResourceNote());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_rsid, m_chrPosition, m_genePosition, m_proteinNote, m_resourceNote, m_isInDel);
  }
}
