package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;

import static org.pharmgkb.pharmcat.reporter.TextConstants.GENOTYPE_DELIMITER;


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

  public boolean isRecommendation() {
    return m_recommendation;
  }
}
