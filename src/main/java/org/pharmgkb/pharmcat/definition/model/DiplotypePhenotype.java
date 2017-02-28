package org.pharmgkb.pharmcat.definition.model;

import java.util.List;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * Model that represents a functional diplotype (e.g. No function/No function) and a related phenotype
 * (e.g. Poor Metabolizer)
 *
 * @author Ryan Whaley
 */
public class DiplotypePhenotype {

  @SerializedName("phenotype")
  @Expose
  private String m_phenotype;
  @SerializedName("diplotype")
  @Expose
  private List<String> m_diplotype;

  public String getPhenotype() {
    return m_phenotype;
  }

  public void setPhenotype(String phenotype) {
    m_phenotype = phenotype;
  }

  public List<String> getDiplotype() {
    return m_diplotype;
  }

  public void setDiplotype(List<String> diplotype) {
    m_diplotype = diplotype;
  }
}
