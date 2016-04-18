package org.pharmgkb.pharmcat.haplotype.model.json;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * Root JSON result object.
 */
public class HaplotyperResult {
  @SerializedName("metadata")
  @Expose
  private Metadata m_metadata;
  @SerializedName("diplotypeCalls")
  @Expose
  private List<DiplotypeCall> m_diplotypeCalls = new ArrayList<>();


  public Metadata getMetadata() {
    return m_metadata;
  }

  public void setMetadata(Metadata metadata) {
    m_metadata = metadata;
  }


  public List<DiplotypeCall> getDiplotypeCalls() {
    return m_diplotypeCalls;
  }

  public void addDiplotypeCall(DiplotypeCall call) {
    m_diplotypeCalls.add(call);
  }
}
