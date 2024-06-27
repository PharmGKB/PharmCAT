package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * PharmGKB Guideline Annotation Model.
 */
public class DosingGuideline {
  @Expose
  @SerializedName("id")
  private String id;
  @Expose
  @SerializedName("name")
  private String name;
  @Expose
  @SerializedName("objCls")
  private String objCls;
  @Expose
  @SerializedName("source")
  private String source;
  @Expose
  @SerializedName("version")
  private Integer m_version;
  @Expose
  @SerializedName("url")
  private String url;
  @Expose
  @SerializedName("relatedChemicals")
  private List<AccessionObject> relatedChemicals = new ArrayList<>();
  @Expose
  @SerializedName("relatedGenes")
  private List<AccessionObject> relatedGenes = new ArrayList<>();
  @Expose
  @SerializedName("relatedAlleles")
  private List<AccessionObject> relatedAlleles = new ArrayList<>();
  @Expose
  @SerializedName("recommendation")
  private boolean m_recommendation;


  public String getId() {
    return id;
  }

  public String getName() {
    return name;
  }

  public String getSource() {
    return source;
  }

  public String getObjCls() {
    return objCls;
  }

  public Integer getVersion() {
    return m_version;
  }

  public String getUrl() {
    return url;
  }

  public List<AccessionObject> getRelatedChemicals() {
    return relatedChemicals;
  }

  public List<AccessionObject> getRelatedGenes() {
    return relatedGenes;
  }

  public List<AccessionObject> getRelatedAlleles() {
    return relatedAlleles;
  }

  public boolean isRecommendation() {
    return m_recommendation;
  }
}
