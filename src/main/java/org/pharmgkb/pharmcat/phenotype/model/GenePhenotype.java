package org.pharmgkb.pharmcat.phenotype.model;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.SortedSet;
import java.util.stream.Collectors;
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
  private SortedSet<DiplotypeRecord> m_diplotypes;
  @SerializedName("namedAlleles")
  @Expose
  private List<HaplotypeRecord> m_namedAlleles;
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


  /**
   * Map of a haplotype name (e.g. *1) to an activity value (e.g. 1.0) if applicable. Not all genes use activity values
   */
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
    return !m_activityValues.isEmpty() && m_activityValues.values().stream()
        .filter(StringUtils::isNotBlank)
        .anyMatch(v -> !v.equalsIgnoreCase(TextConstants.NA));
  }


  public String getHaplotypeFunction(String haplotype) {
    if (StringUtils.isBlank(haplotype) || m_namedAlleles == null) {
      return UNASSIGNED_FUNCTION;
    }

    return m_namedAlleles.stream()
        .filter(n -> n.getName().equals(haplotype))
        .findFirst()
        .map(HaplotypeRecord::getFunctionValue)
        .orElse(UNASSIGNED_FUNCTION);
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
  public SortedSet<DiplotypeRecord> getDiplotypes() {
    return m_diplotypes;
  }

  public List<HaplotypeRecord> getNamedAlleles() {
    return m_namedAlleles;
  }

  public Optional<DiplotypeRecord> findDiplotype(Map<String,Integer> diplotypeKey) {
    List<DiplotypeRecord> diplotypes = m_diplotypes.stream()
        .filter(d -> d.matchesKey(diplotypeKey))
        .toList();
    if (diplotypes.size() == 1) {
      return Optional.of(diplotypes.get(0));
    } else if (diplotypes.isEmpty()) {
      return Optional.empty();
    }
    // should never happen, DataManager should have caught this
    throw new IllegalStateException(diplotypes.size() + " diplotypes found for " + m_gene + " for " +
        diplotypeKey.keySet().stream()
            .map(k -> k + " (" + diplotypeKey.get(k) + ")")
            .collect(Collectors.joining(", ")));
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
