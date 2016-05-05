package org.pharmgkb.pharmcat.definition.model;

import java.util.Objects;
import com.google.gson.annotations.SerializedName;


/**
 * All the information used to describe a particular location of a variant.
 *
 * @author Ryan Whaley
 */
public class VariantLocus {
  @SerializedName("position")
  private int m_position;
  @SerializedName("rsid")
  private String m_rsid;
  @SerializedName("chromosomeHgvsName")
  private String m_chromosomeHgvsName;
  @SerializedName("geneHgvsName")
  private String m_geneHgvsName;
  @SerializedName("proteinNote")
  private String m_proteinNote;
  @SerializedName("resourceNote")
  private String m_resourceNote;
  @SerializedName("isInDel")
  private boolean m_isInDel;


  public VariantLocus(int position, String chromosomeHgvsName) {
    m_position = position;
    m_chromosomeHgvsName = chromosomeHgvsName;
  }


  /**
   * The (start) position on the chromosomal sequence.
   */
  public int getPosition() {
    return m_position;
  }

  /**
   * The name use for this location on the chromosomal sequence, should be relative to plus strand
   */
  public String getChromosomeHgvsName() {
    return m_chromosomeHgvsName;
  }


  /**
   * The name use for this location on the gene sequence, relative to the strand the gene is on
   */
  public String getGeneHgvsName() {
    return m_geneHgvsName;
  }

  public void setGeneHgvsName(String geneHgvsName) {
    m_geneHgvsName = geneHgvsName;
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


  /**
   * A common name for this variant, usually specified by some specialized resource
   */
  public String getResourceNote() {
    return m_resourceNote;
  }

  public void setResourceNote(String resourceNote) {
    m_resourceNote = resourceNote;
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
        m_position == that.getPosition() &&
        Objects.equals(m_chromosomeHgvsName, that.getChromosomeHgvsName()) &&
        Objects.equals(m_geneHgvsName, that.getGeneHgvsName()) &&
        Objects.equals(m_proteinNote, that.getProteinNote()) &&
        Objects.equals(m_rsid, that.getRsid()) &&
        Objects.equals(m_resourceNote, that.getResourceNote());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_position, m_chromosomeHgvsName, m_geneHgvsName, m_proteinNote, m_rsid, m_resourceNote, m_isInDel);
  }
}
