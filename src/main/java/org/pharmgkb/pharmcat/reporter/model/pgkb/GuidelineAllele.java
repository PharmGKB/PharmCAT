package org.pharmgkb.pharmcat.reporter.model.pgkb;

import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * PharmGKB Guideline Annotation Allele Model
 */
@Deprecated
public class GuidelineAllele {
  @SerializedName("id")
  @Expose
  private Integer m_id;
  @SerializedName("_label")
  @Expose
  private String m_label;
  @SerializedName("function")
  @Expose
  private OntologyTerm m_functionTerm;

  public Integer getId() {
    return m_id;
  }

  public void setId(Integer id) {
    m_id = id;
  }

  public String getLabel() {
    return m_label;
  }

  public void setLabel(String label) {
    m_label = label;
  }

  public OntologyTerm getFunctionTerm() {
    return m_functionTerm;
  }

  public void setFunctionTerm(OntologyTerm functionTerm) {
    m_functionTerm = functionTerm;
  }
}
