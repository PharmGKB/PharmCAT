package org.pharmgkb.pharmcat.definition.model;

/**
 * All the information used to describe a particular location of a Variant
 *
 * @author Ryan Whaley
 */
public class VariantLocus {

  private String m_resourceNote;
  private String m_proteinNote;
  private String m_chrPosition;
  private String m_genePosition;
  private String m_rsid;

  public String toString() {
    return m_chrPosition;
  }

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
   * The name use for this location on the protein sequence
   */
  public String getProteinNote() {
    return m_proteinNote;
  }

  public void setProteinNote(String proteinNote) {
    m_proteinNote = proteinNote;
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
   * The identifier use for this location from dbSNP
   */
  public String getRsid() {
    return m_rsid;
  }

  public void setRsid(String rsid) {
    m_rsid = rsid;
  }
}
