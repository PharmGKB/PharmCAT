package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class CPICinteraction {

  @SerializedName("objCls")
  @Expose
  private String objCls;
  @SerializedName("@id")
  @Expose
  private String Id;
  @SerializedName("@context")
  @Expose
  private String Context;
  @SerializedName("id")
  @Expose
  private String id;
  @SerializedName("name")
  @Expose
  private String name;
  @SerializedName("groups")
  @Expose
  private List<Group> groups = new ArrayList<Group>();
  @SerializedName("relatedChemicals")
  @Expose
  private List<RelatedChemical> relatedChemicals = new ArrayList<RelatedChemical>();
  @SerializedName("relatedGenes")
  @Expose
  private List<RelatedGene> relatedGenes = new ArrayList<RelatedGene>();
  @SerializedName("source")
  @Expose
  private String source;
  @SerializedName("summaryHtml")
  @Expose
  private String summaryHtml;
  @SerializedName("textHtml")
  @Expose
  private String textHtml;

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
   *     The Id
   *
  public String getId() {
  return Id;
  }

  /**
   *
   * @param Id
   *     The @id
   *
  public void setId(String Id) {
  this.Id = Id;
  }

  /**
   *
   * @return
   *     The Context
   *
  public String getContext() {
  return Context;
  }

  /**
   *
   * @param Context
   *     The @context
   */
  public void setContext(String Context) {
    this.Context = Context;
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
   *     The groups
   */
  public List<Group> getGroups() {
    return groups;
  }

  /**
   *
   * @param groups
   *     The groups
   */
  public void setGroups(List<Group> groups) {
    this.groups = groups;
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
    return summaryHtml;
  }

  /**
   *
   * @param summaryHtml
   *     The summaryHtml
   */
  public void setSummaryHtml(String summaryHtml) {
    this.summaryHtml = summaryHtml;
  }

  /**
   *
   * @return
   *     The textHtml
   */
  public String getTextHtml() {
    return textHtml;
  }

  /**
   *
   * @param textHtml
   *     The textHtml
   */
  public void setTextHtml(String textHtml) {
    this.textHtml = textHtml;
  }
}
