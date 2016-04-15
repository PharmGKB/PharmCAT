
package org.cpic.haplotype.model.json;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class Variant {
  @SerializedName("m_position")
  @Expose
  private int m_position;
  @SerializedName("m_rsid")
  @Expose
  private String m_rsid;
  @SerializedName("m_vcfCall")
  @Expose
  private String m_vcfCall;


  public Variant(int pos, String rsids, String call) {
    m_position = pos;
    m_rsid = rsids;
    m_vcfCall = call;
  }


  public int getPosition() {
    return m_position;
  }

  public String getRsid() {
    return m_rsid;
  }

  public String getVcfCall() {
    return m_vcfCall;
  }
}
