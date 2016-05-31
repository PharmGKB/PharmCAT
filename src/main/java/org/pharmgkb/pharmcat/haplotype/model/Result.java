package org.pharmgkb.pharmcat.haplotype.model;

import java.util.ArrayList;
import java.util.List;
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
  private List<GeneCall> m_geneCalls = new ArrayList<>();


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
}
