
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.api.client.repackaged.com.google.common.base.Joiner;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.SampleAllele;


public class Variant implements Comparable<Variant>  {
  private static final String sf_rsidFormat = "%s%s";
  private static final String sf_positionFormat = "%d%s";
  private static final Joiner sf_vcfAlleleJoiner = Joiner.on(",");

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


  public Variant(@Nonnull VariantLocus variant, @Nonnull SampleAllele allele) {
    String call;
    String vcfAlleles = sf_vcfAlleleJoiner.join(allele.getVcfAlleles());
    if (allele.isPhased()) {
      call = allele.getAllele1() + "|" + allele.getAllele2();
    } else {
      call = allele.getAllele1() + "/" + allele.getAllele2();
    }
    initialize(variant.getPosition(), variant.getRsid(), call, variant.getVcfPosition(), vcfAlleles);
  }


  public Variant(int pos, @Nullable String rsids, @Nullable String call, int vcfPosition, @Nullable String vcfAlleles) {
    initialize(pos, rsids, call, vcfPosition, vcfAlleles);
  }

  private void initialize(int pos, @Nullable String rsids, @Nullable String call, int vcfPosition, @Nullable String vcfAlleles) {
    m_position = pos;
    m_rsid = rsids;
    m_vcfCall = call;
    if (call != null) {
      m_isPhased = call.contains("|");
    }
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
    return Comparator.nullsLast(String::compareTo).compare(m_vcfCall, o.getVcfCall());
  }

  public String printDisplay() {
    String[] alleles = getVcfCall().split("[|\\/]");
    if (m_rsid != null) {
      return Arrays.stream(alleles).map(a -> getRsid()+a).collect(Collectors.joining("/"));
    }
    else {
      return Arrays.stream(alleles).map(a -> getPosition()+a).collect(Collectors.joining("/"));
    }
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
