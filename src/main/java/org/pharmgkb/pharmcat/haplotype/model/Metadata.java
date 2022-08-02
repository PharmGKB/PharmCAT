
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Date;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


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
  private Date m_timetamp;
  @Expose
  @SerializedName("topCandidatesOnly")
  private boolean m_topCandidatesOnly;
  @Expose
  @SerializedName("findCombinations")
  private boolean m_findCombinations;
  @Expose
  @SerializedName("callCyp2d")
  private boolean m_callCyp2d6;


  public Metadata(String namedAlleleMatcherVersion, String genomeBuild, String vcfFilename, Date date,
      boolean topCandidatesOnly, boolean findCombinations, boolean callCyp2d6) {
    m_namedAlleleMatcherVersion = namedAlleleMatcherVersion;
    m_genomeBuild = genomeBuild;
    m_inputFilename = vcfFilename;
    m_timetamp = date;
    m_topCandidatesOnly = topCandidatesOnly;
    m_findCombinations = findCombinations;
    m_callCyp2d6 = callCyp2d6;
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

  public Date getTimetamp() {
    return m_timetamp;
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
}
