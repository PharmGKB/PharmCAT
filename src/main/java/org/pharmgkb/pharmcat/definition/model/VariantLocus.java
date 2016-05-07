package org.pharmgkb.pharmcat.definition.model;

import java.util.Objects;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * All the information used to describe a particular location of a variant.
 *
 * @author Ryan Whaley
 */
public class VariantLocus implements Comparable<VariantLocus> {
  @Expose
  @SerializedName("position")
  private int m_position;
  @Expose
  @SerializedName("rsid")
  private String m_rsid;
  @Expose
  @SerializedName("chromosomeHgvsName")
  private String m_chromosomeHgvsName;
  @Expose
  @SerializedName("geneHgvsName")
  private String m_geneHgvsName;
  @Expose
  @SerializedName("proteinNote")
  private String m_proteinNote;
  @Expose
  @SerializedName("resourceNote")
  private String m_resourceNote;
  @Expose
  @SerializedName("isInDel")
  private boolean m_isInDel;
  @Expose
  @SerializedName("isRepeat")
  private boolean m_isRepeat;


  public VariantLocus(int position, @Nonnull String chromosomeHgvsName) {
    Preconditions.checkNotNull(chromosomeHgvsName);
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
  public @Nonnull String getChromosomeHgvsName() {
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


  /**
   * Gets whether this locus has an insertion or deletion.
   */
  public boolean isInDel() {
    return m_isInDel;
  }

  public void setInDel(boolean inDel) {
    Preconditions.checkState(!m_isInDel, "Cannot be both indel and repeat");
    m_isInDel = inDel;
  }

  /**
   * Gets whether this locus has a repeat.
   */
  public boolean isRepeat() {
    return m_isRepeat;
  }

  public void setRepeat(boolean repeat) {
    Preconditions.checkState(!m_isRepeat, "Cannot be both indel and repeat");
    m_isRepeat = repeat;
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
    return m_position == that.getPosition() &&
        m_isInDel == that.isInDel() &&
        m_isRepeat == that.isRepeat() &&
        Objects.equals(m_chromosomeHgvsName, that.getChromosomeHgvsName()) &&
        Objects.equals(m_geneHgvsName, that.getGeneHgvsName()) &&
        Objects.equals(m_proteinNote, that.getProteinNote()) &&
        Objects.equals(m_rsid, that.getRsid()) &&
        Objects.equals(m_resourceNote, that.getResourceNote());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_position, m_chromosomeHgvsName, m_geneHgvsName, m_proteinNote, m_rsid, m_resourceNote);
  }


  @Override
  public int compareTo(@Nonnull VariantLocus o) {

    int rez = Integer.compare(m_position, o.getPosition());
    if (rez != 0) {
      return rez;
    }
    return m_chromosomeHgvsName.compareTo(o.getChromosomeHgvsName());
  }
}
