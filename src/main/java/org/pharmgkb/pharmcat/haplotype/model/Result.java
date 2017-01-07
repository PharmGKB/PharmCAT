package org.pharmgkb.pharmcat.haplotype.model;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.haplotype.SampleAllele;


/**
 * Root JSON result object.
 */
public class Result {
  @SerializedName("metadata")
  @Expose
  private Metadata m_metadata;

  @SerializedName("results")
  @Expose
  private List<GeneCall> m_geneCalls = new ArrayList<>();

  @Expose
  @SerializedName("annotatedAlleles")
  private Set<SampleAllele> m_annotatedAlleles;


  public Metadata getMetadata() {
    return m_metadata;
  }

  public void setMetadata(Metadata metadata) {
    m_metadata = metadata;
  }


  public List<GeneCall> getGeneCalls() {
    return m_geneCalls;
  }

  public void addDiplotypeCall(GeneCall call) {
    m_geneCalls.add(call);
  }


  public Set<SampleAllele> getAnnotatedAlleles() {
    return m_annotatedAlleles;
  }

  public void setAnnotatedAlleles(Set<SampleAllele> alleles) {
    m_annotatedAlleles = alleles;
  }
}
