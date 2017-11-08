package org.pharmgkb.pharmcat.reporter.model;

import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * Message annotations to put on specific variants
 *
 * @author Ryan Whaley
 */
public class MessageVariant {

  private static final int sf_rowLength = 8;

  @Expose
  @SerializedName("rsid")
  private String m_rsid;
  @Expose
  @SerializedName("chr")
  private String m_chr;
  @Expose
  @SerializedName("position")
  private String m_position;
  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("guideline")
  private String m_guideline;

  public MessageVariant(@Nonnull String row) throws RuntimeException {
    String[] fields = row.split("\\t");
    Preconditions.checkArgument(fields.length == sf_rowLength,
        "Message Variant column count mismatch, got " + fields.length + ", expected " + sf_rowLength);

    setGene(fields[1]);
    setGuideline(fields[2]);
    setRsid(fields[3]);
    setChr(fields[4]);
    setPosition(fields[5]);
  }

  public String getRsid() {
    return m_rsid;
  }

  public void setRsid(String rsid) {
    m_rsid = rsid;
  }

  public String getChr() {
    return m_chr;
  }

  public void setChr(String chr) {
    m_chr = chr;
  }

  public String getPosition() {
    return m_position;
  }

  public void setPosition(String position) {
    m_position = position;
  }

  public String getGene() {
    return m_gene;
  }

  public void setGene(String gene) {
    m_gene = gene;
  }

  public String getGuideline() {
    return m_guideline;
  }

  public void setGuideline(String guideline) {
    m_guideline = guideline;
  }

  public String toString() {
    if (m_rsid != null) {
      return m_rsid;
    }
    return getChr() + ":" + m_position;
  }
}
