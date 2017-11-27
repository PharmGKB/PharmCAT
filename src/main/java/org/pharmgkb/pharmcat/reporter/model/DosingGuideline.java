package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.result.Markdown;


public class DosingGuideline {

  @SerializedName("objCls")
  @Expose
  private String objCls;
  @SerializedName("id")
  @Expose
  private String id;
  @SerializedName("name")
  @Expose
  private String name;
  @SerializedName("relatedChemicals")
  @Expose
  private List<RelatedChemical> relatedChemicals = new ArrayList<>();
  @SerializedName("relatedGenes")
  @Expose
  private List<RelatedGene> relatedGenes = new ArrayList<>();
  @SerializedName("source")
  @Expose
  private String source;
  @SerializedName("summaryMarkdown")
  @Expose
  private Markdown m_summaryMarkdown;
  @SerializedName("textMarkdown")
  @Expose
  private Markdown m_textMarkdown;
  @SerializedName("history")
  @Expose
  private List<History> m_history;
  @SerializedName("recommendation")
  @Expose
  private boolean m_recommendation;

  /**
   *
   * @return
   *     The objCls
   */
  public String getObjCls() {
    return objCls;
  }

  /**
   *
   * @param objCls
   *     The objCls
   */
  public void setObjCls(String objCls) {
    this.objCls = objCls;
  }

  /**
   *
   * @return
   *     The id
   */
  public String getId() {
    return id;
  }

  /**
   *
   * @param id
   *     The id
   */
  public void setId(String id) {
    this.id = id;
  }

  /**
   *
   * @return
   *     The name
   */
  public String getName() {
    return name;
  }

  /**
   *
   * @param name
   *     The name
   */
  public void setName(String name) {
    this.name = name;
  }

  /**
   *
   * @return
   *     The relatedChemicals
   */
  public List<RelatedChemical> getRelatedChemicals() {
    return relatedChemicals;
  }

  /**
   *
   * @param relatedChemicals
   *     The relatedChemicals
   */
  public void setRelatedChemicals(List<RelatedChemical> relatedChemicals) {
    this.relatedChemicals = relatedChemicals;
  }

  /**
   *
   * @return
   *     The relatedGenes
   */
  public List<RelatedGene> getRelatedGenes() {
    return relatedGenes;
  }

  /**
   *
   * @param relatedGenes
   *     The relatedGenes
   */
  public void setRelatedGenes(List<RelatedGene> relatedGenes) {
    this.relatedGenes = relatedGenes;
  }

  /**
   *
   * @return
   *     The source
   */
  public String getSource() {
    return source;
  }

  /**
   *
   * @param source
   *     The source
   */
  public void setSource(String source) {
    this.source = source;
  }

  /**
   *
   * @return
   *     The summaryHtml
   */
  public String getSummaryHtml() {
    return m_summaryMarkdown.getHtml();
  }

  /**
   *
   * @return
   *     The textHtml
   */
  public String getTextHtml() {
    return m_textMarkdown.getHtml();
  }

  public Markdown getSummaryMarkdown() {
    return m_summaryMarkdown;
  }

  public void setSummaryMarkdown(Markdown summaryMarkdown) {
    m_summaryMarkdown = summaryMarkdown;
  }

  public Markdown getTextMarkdown() {
    return m_textMarkdown;
  }

  public void setTextMarkdown(Markdown textMarkdown) {
    m_textMarkdown = textMarkdown;
  }

  public List<History> getHistory() {
    return m_history;
  }

  public boolean isRecommendation() {
    return m_recommendation;
  }

  public void setRecommendation(boolean recommendation) {
    m_recommendation = recommendation;
  }
}
