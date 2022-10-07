package org.pharmgkb.pharmcat.phenotype.model;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.ActivityUtils;


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

  public void setGene(String gene) {
    m_gene = gene;
  }

  /**
   * Map of haplotype name (e.g. *1) to a phenotype (e.g. Increased function)
   */
  public Map<String, String> getHaplotypes() {
    return m_haplotypes;
  }

  public void setHaplotypes(Map<String, String> haplotypes) {
    m_haplotypes = haplotypes;
  }


  public Map<String,String> getActivityValues() {
    return m_activityValues;
  }

  public void setActivityValues(Map<String, String> activityValues) {
    m_activityValues = activityValues;
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

  public void setDiplotypes(List<DiplotypeRecord> diplotypes) {
    m_diplotypes = diplotypes;
  }

  /**
   * *1/*3 to phenotype
   *
   * @param diplotype a String like "*1/*4"
   * @return a phenotype value like "Normal Metabolizer"
   */
  public String getPhenotypeForDiplotype(Diplotype diplotype) {
    if (diplotype.isUnknown()) {
      return TextConstants.NO_RESULT;
    }

    Set<String> phenos = getDiplotypes().stream()
        .filter(d -> d.getDiplotypeKey().equals(diplotype.makeLookupMap()))
        .map(DiplotypeRecord::getGeneResult)
        .collect(Collectors.toSet());
    if (phenos.size()>1) {
      throw new IllegalStateException("More than one phenotype match made for " + getGene() + " " + diplotype + ": " + String.join("; ", phenos));
    } else if (phenos.size() == 0) {
      if (diplotype.getAllele2() == null && diplotype.getGene().equals("DPYD")) {
        return "";
      }
      return TextConstants.NA;
    } else {
      return phenos.iterator().next();
    }
  }

  /**
   * Get the phenotype value for the activity score of a given diplotype.
   * @param diplotype a diplotype to analyze, should have an activity score
   * @return the phenotype that matches the phenotype of the diplotype
   */
  public String getPhenotypeForActivity(Diplotype diplotype) {
    if (!isMatchedByActivityScore()) {
      return null;
    }
    return getDiplotypes().stream()
        .filter(d -> d.getLookupKey().equals(diplotype.getActivityScore()))
        .findFirst()
        .map(DiplotypeRecord::getGeneResult)
        .orElse(TextConstants.NO_RESULT);
  }

  /**
   * Gets the activity score for the given diplotype using the alleles of the diplotype.
   *
   * @param diplotype a diplotype which should have alleles
   * @return the activity score string for the given diplotype
   */
  public String getActivityForDiplotype(Diplotype diplotype) {
    if (diplotype.isUnknown() || !isMatchedByActivityScore()) {
      return null;
    }

    Set<String> scores = getDiplotypes().stream()
        .filter(d -> d.getDiplotypeKey().equals(diplotype.makeLookupMap()))
        .map(DiplotypeRecord::getLookupKey)
        .collect(Collectors.toSet());
    if (scores.size() > 1) {
      throw new IllegalStateException("More than one phenotype match made for " + getGene() + " " + diplotype + ": " + String.join("; ", scores));
    } else if (scores.size() == 0) {
      return TextConstants.NA;
    } else {
      return ActivityUtils.normalize(scores.iterator().next());
    }
  }

  /**
   * Gets the lookup key for the given bare diplotype of this gene.
   *
   * Runs through all the possible diplotypes in this {@link GenePhenotype} record to find one that matches
   *
   * @param diplotype in the form of "*1/*3"
   * @return the lookup key related to this diplotype
   */
  public String getLookupKeyForDiplotype(Diplotype diplotype) {
    if (diplotype.isUnknown()) {
      return TextConstants.NO_RESULT;
    }

    Set<String> keys = m_diplotypes.stream()
        .filter(d -> d.getDiplotypeKey().equals(diplotype.makeLookupMap()))
        .map(DiplotypeRecord::getLookupKey)
        .collect(Collectors.toSet());
    if (keys.size() > 1) {
      throw new IllegalStateException("More than one key match made for " + getGene() + " " + diplotype + ": " + String.join("; ", keys));
    } else if (keys.size() == 0) {
      return TextConstants.NA;
    } else {
      return keys.iterator().next();
    }
  }


  /**
   * Gets the version of the {@link DataSource} the phenotype mapping is from.
   */
  public String getVersion() {
    return m_version;
  }

  public void setVersion(String version) {
    m_version = version;
  }


  @Override
  public String toString() {
    return m_gene;
  }
}
