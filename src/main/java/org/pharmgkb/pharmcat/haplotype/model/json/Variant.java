
package org.pharmgkb.pharmcat.haplotype.model.json;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class Variant implements Comparable<Variant>  {
  @SerializedName("position")
  @Expose
  private int m_position;
  @SerializedName("rsid")
  @Expose
  private String m_rsid;
  @SerializedName("vcfCall")
  @Expose
  private String m_vcfCall;


  public Variant(int pos, @Nullable String rsids, @Nonnull String call) {
    Preconditions.checkNotNull(call);
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

  @Override
  public int compareTo(@Nonnull Variant o) {

    int rez = Integer.compare(m_position, o.getPosition());
    if (rez != 0) {
      return rez;
    }
    return m_vcfCall.compareTo(o.getVcfCall());
  }
}
