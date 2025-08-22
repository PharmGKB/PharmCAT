
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Date;
import java.util.Map;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.jspecify.annotations.Nullable;


public class Metadata {
  @Expose
  @SerializedName("namedAlleleMatcherVersion")
  private String m_namedAlleleMatcherVersion;
  @Expose
  @SerializedName("genomeBuild")
  private String m_genomeBuild;
  @Expose
  @SerializedName("inputFilename")
  private String m_inputFilename;
  @Expose
  @SerializedName("timestamp")
  private Date m_timestamp;
  @Expose
  @SerializedName("topCandidatesOnly")
  private boolean m_topCandidatesOnly;
  @Expose
  @SerializedName("findCombinations")
  private boolean m_findCombinations;
  @Expose
  @SerializedName("callCyp2d")
  private boolean m_callCyp2d6;
  @Expose
  @SerializedName("sampleId")
  private String m_sampleId;
  @Expose
  @SerializedName("sampleProps")
  private Map<String, String> m_sampleProps;


  public Metadata(String namedAlleleMatcherVersion, String genomeBuild, String vcfFilename, Date date,
      boolean topCandidatesOnly, boolean findCombinations, boolean callCyp2d6, String sampleId) {
    m_namedAlleleMatcherVersion = namedAlleleMatcherVersion;
    m_genomeBuild = genomeBuild;
    m_inputFilename = vcfFilename;
    m_timestamp = date;
    m_topCandidatesOnly = topCandidatesOnly;
    m_findCombinations = findCombinations;
    m_callCyp2d6 = callCyp2d6;
    m_sampleId = sampleId;
  }


  public String getNamedAlleleMatcherVersion() {
    return m_namedAlleleMatcherVersion;
  }

  public String getGenomeBuild() {
    return m_genomeBuild;
  }

  public String getInputFilename() {
    return m_inputFilename;
  }

  public Date getTimestamp() {
    return m_timestamp;
  }

  public boolean isTopCandidatesOnly() {
    return m_topCandidatesOnly;
  }

  public boolean isFindCombinations() {
    return m_findCombinations;
  }

  public boolean isCallCyp2d6() {
    return m_callCyp2d6;
  }


  public String getSampleId() {
    return m_sampleId;
  }

  public void setSampleId(String sampleId) {
    m_sampleId = sampleId;
  }

  public @Nullable Map<String, String> getSampleProps() {
    return m_sampleProps;
  }

  public void setSampleProps(Map<String, String> sampleProps) {
    m_sampleProps = sampleProps;
  }
}
