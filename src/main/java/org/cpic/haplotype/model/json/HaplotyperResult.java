package org.cpic.haplotype.model.json;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class HaplotyperResult {
  @SerializedName("diplotypeCalls")
  @Expose
  private List<DiplotypeCall> m_diplotypeCalls = new ArrayList<>();
  @SerializedName("metadata")
  @Expose
  private Metadata m_metadata;


  public List<DiplotypeCall> getDiplotypeCalls() {
    return m_diplotypeCalls;
  }

  public void addDiplotypeCall(DiplotypeCall call) {
    m_diplotypeCalls.add(call);
  }


  public Metadata getMetadata() {
    return m_metadata;
  }

  public void setMetadata(Metadata metadata) {
    this.m_metadata = metadata;
  }
}
