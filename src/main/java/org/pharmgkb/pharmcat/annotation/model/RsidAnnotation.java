package org.pharmgkb.pharmcat.annotation.model;

import javax.annotation.Nonnull;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * A reporting annotation for an RSID.
 *
 * @author Mark Woon
 */
public class RsidAnnotation {
  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("rsid")
  private String m_rsid;


  public RsidAnnotation(@Nonnull String gene, @Nonnull String rsid) {
    m_gene = gene;
    m_rsid = rsid;
  }


  public String getGene() {
    return m_gene;
  }

  public String getRsid() {
    return m_rsid;
  }
}
