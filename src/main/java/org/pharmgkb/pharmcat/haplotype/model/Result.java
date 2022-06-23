package org.pharmgkb.pharmcat.haplotype.model;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * Root JSON result object.
 */
public class Result {
  @SerializedName("metadata")
  @Expose
  private Metadata m_metadata;

  @SerializedName("results")
  @Expose
  private final List<GeneCall> m_geneCalls = new ArrayList<>();

  @SerializedName("vcfWarnings")
  @Expose
  private Map<String, Collection<String>> m_vcfWarnings;


  public Metadata getMetadata() {
    return m_metadata;
  }

  public void setMetadata(Metadata metadata) {
    m_metadata = metadata;
  }


  public List<GeneCall> getGeneCalls() {
    return m_geneCalls;
  }

  public void addGeneCall(GeneCall call) {
    m_geneCalls.add(call);
  }


  /**
   * Gets warnings from reading VCF data, keyed to chromosomal position.
   */
  public Map<String, Collection<String>> getVcfWarnings() {
    return m_vcfWarnings;
  }

  public void setVcfWarnings(Map<String, Collection<String>> warnings) {
    m_vcfWarnings = warnings;
  }
}
