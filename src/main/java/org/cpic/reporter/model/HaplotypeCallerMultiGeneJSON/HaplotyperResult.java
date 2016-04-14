package org.cpic.reporter.model.HaplotypeCallerMultiGeneJSON;

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
  private Foo m_metadata;


  public List<DiplotypeCall> getDiplotypeCalls() {
    return m_diplotypeCalls;
  }

  public void addDiplotypeCall(DiplotypeCall call) {
    m_diplotypeCalls.add(call);
  }


  public Foo getMetadata() {
    return m_metadata;
  }

  public void setMetadata(Foo metadata) {
    this.m_metadata = metadata;
  }
}
