package org.pharmgkb.pharmcat.definition.model;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;


/**
 * Model object that contains all the phenotype mapping data needed for a gene
 *
 * @author Ryan Whaley
 */
public class GenePhenotype {
  private static final String NO_RESULT = "No Result";

  @SerializedName("gene")
  @Expose
  private String m_gene;
  @SerializedName("haplotypes")
  @Expose
  private Map<String,String> m_haplotypes;
  @SerializedName("diplotypes")
  @Expose
  private List<DiplotypeRecord> m_diplotypes;

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


  @Nullable
  public String lookupHaplotype(@Nullable String hap) {
    if (hap == null) {
      return null;
    }
    return m_haplotypes.get(hap);
  }

  public Optional<String> findHaplotypeFunction(String haplotype) {
    if (StringUtils.isBlank(haplotype)) {
      return Optional.empty();
    }
    return Optional.ofNullable(m_haplotypes.get(haplotype));
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
   * @return a phenotype value Normal Metabolizer
   */
  public String getPhenotypeForDiplotype(Diplotype diplotype) {
    if (diplotype.isUnknown()) {
      return NO_RESULT;
    }

    Set<String> phenos = getDiplotypes().stream()
        .filter(d -> d.getDiplotypeKey().equals(diplotype.makeLookupMap()))
        .map(DiplotypeRecord::getGeneresult)
        .collect(Collectors.toSet());
    if (phenos.size()>1) {
      throw new IllegalStateException("More than one phenotype match made for " + getGene() + " " + diplotype + ": " + String.join("; ", phenos));
    } else if (phenos.size() == 0) {
      return "N/A";
    } else {
      return phenos.iterator().next();
    }
  }

  /**
   * Gets the lookup key for the given bare diplotype of this gene.
   * @param diplotype in the form of "*1/*3"
   * @return the lookup key related to this diplotype
   */
  public String getLookupKeyForDiplotype(Diplotype diplotype) {
    if (diplotype.isUnknown()) {
      return NO_RESULT;
    }

    Set<String> keys = m_diplotypes.stream()
        .filter(d -> d.getDiplotypeKey().equals(diplotype.makeLookupMap()))
        .map(DiplotypeRecord::getLookupKey)
        .collect(Collectors.toSet());
    if (keys.size() > 1) {
      throw new IllegalStateException("More than one key match made for " + getGene() + " " + diplotype + ": " + String.join("; ", keys));
    } else if (keys.size() == 0) {
      return "N/A";
    } else {
      return keys.iterator().next();
    }
  }

  public String toString() {
    return m_gene;
  }
}
