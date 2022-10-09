package org.pharmgkb.pharmcat.phenotype.model;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


/**
 * Model object that contains all the phenotype mapping data needed for a gene.
 *
 * @author Ryan Whaley
 */
public class GenePhenotype {
  public static final String UNASSIGNED_FUNCTION = "Unassigned function";

  @SerializedName("gene")
  @Expose
  private String m_gene;
  @SerializedName("haplotypes")
  @Expose
  private Map<String, String> m_haplotypes;
  @SerializedName("activityValues")
  @Expose
  private Map<String,String> m_activityValues = new HashMap<>();
  @SerializedName("diplotypes")
  @Expose
  private List<DiplotypeRecord> m_diplotypes;
  @SerializedName(value = "version", alternate = {"cpicVersion"})
  @Expose
  private String m_version;


  /**
   * The HGNC gene symbol
   */
  public String getGene() {
    return m_gene;
  }


  /**
   * Map of haplotype name (e.g. *1) to a phenotype (e.g. Increased function)
   */
  public Map<String, String> getHaplotypes() {
    return m_haplotypes;
  }


  public Map<String,String> getActivityValues() {
    return m_activityValues;
  }

  public void assignActivity(Haplotype haplotype) {
    if (haplotype == null || haplotype.isUnknown()) {
      return;
    }
    haplotype.setActivityValue(m_activityValues.getOrDefault(haplotype.getName(), TextConstants.NA));
  }

  public boolean isMatchedByActivityScore() {
    return m_activityValues.size() > 0;
  }


  public String getHaplotypeFunction(String haplotype) {
    if (StringUtils.isBlank(haplotype) || m_haplotypes == null) {
      return UNASSIGNED_FUNCTION;
    }
    return m_haplotypes.getOrDefault(haplotype, UNASSIGNED_FUNCTION);
  }

  public @Nullable String getHaplotypeActivity(String haplotype) {
    return m_activityValues.get(haplotype);
  }

  public @Nullable Float getHaplotypeActivityScore(String haplotype) {
    String value = m_activityValues.get(haplotype);
    if (value == null || value.equalsIgnoreCase(TextConstants.NA)) {
      return null;
    }
    return Float.valueOf(value);
  }


  /**
   * List of all diplotype to phenotype mappings for this gene
   */
  public List<DiplotypeRecord> getDiplotypes() {
    return m_diplotypes;
  }


  /**
   * Gets the version of the {@link DataSource} the phenotype mapping is from.
   */
  public String getVersion() {
    return m_version;
  }


  @Override
  public String toString() {
    return m_gene;
  }
}
