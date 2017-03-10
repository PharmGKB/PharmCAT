
package org.pharmgkb.pharmcat.haplotype.model;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class Variant implements Comparable<Variant>  {
  private static final String sf_rsidFormat = "%s%s";
  private static final String sf_positionFormat = "%d%s";

  @Expose
  @SerializedName("position")
  private int m_position;
  @Expose
  @SerializedName("rsid")
  private String m_rsid;
  @Expose
  @SerializedName("vcfCall")
  private String m_vcfCall;
  private boolean m_isPhased;
  private int m_vcfPosition;
  private String m_vcfAlleles;


  public Variant(int pos, @Nullable String rsids, @Nonnull String call, int vcfPosition, @Nonnull String vcfAlleles) {
    Preconditions.checkNotNull(call);
    m_position = pos;
    m_rsid = rsids;
    m_vcfCall = call;
    m_isPhased = call.contains("|");
    m_vcfPosition = vcfPosition;
    m_vcfAlleles = vcfAlleles;
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

  public boolean isPhased() {
    return m_isPhased;
  }

  public int getVcfPosition() {
    return m_vcfPosition;
  }

  public String getVcfAlleles() {
    return m_vcfAlleles;
  }

  @Override
  public int compareTo(@Nonnull Variant o) {

    int rez = Integer.compare(m_position, o.getPosition());
    if (rez != 0) {
      return rez;
    }
    rez = Integer.compare(m_vcfPosition, o.getVcfPosition());
    if (rez != 0) {
      return rez;
    }
    return m_vcfCall.compareTo(o.getVcfCall());
  }

  public String toString() {
    if (m_rsid != null) {
      return String.format(sf_rsidFormat, getRsid(), getVcfCall().replaceAll("[|\\/]", ""));
    }
    else {
      return String.format(sf_positionFormat, getPosition(), getVcfCall().replaceAll("[|\\/]", ""));
    }
  }
}
