
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Date;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class Metadata {
  @Expose
  @SerializedName("haplotyperVersion")
  private String m_haplotyperVersion;
  @Expose
  @SerializedName("genomeBuild")
  private String m_genomeBuild;
  @Expose
  @SerializedName("inputFilename")
  private String m_inputFilename;
  @Expose
  @SerializedName("timestamp")
  private Date m_timetamp;


  public Metadata(String haplotyperVersion, String genomeBuild, String vcfFilename, Date date) {
    m_haplotyperVersion = haplotyperVersion;
    m_genomeBuild = genomeBuild;
    m_inputFilename = vcfFilename;
    m_timetamp = date;
  }


  public String getHaplotyperVersion() {
    return m_haplotyperVersion;
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
}
