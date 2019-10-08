package org.pharmgkb.pharmcat.definition.model;

import java.lang.invoke.MethodHandles;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import javax.annotation.Nullable;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.UnexpectedStateException;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Model object that contains all the phenotype mapping data needed for a gene
 *
 * @author Ryan Whaley
 */
public class GenePhenotype {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  @SerializedName("gene")
  @Expose
  private String m_gene;
  @SerializedName("haplotypes")
  @Expose
  private Map<String,String> m_haplotypes;
  @SerializedName("diplotypes")
  @Expose
  private List<DiplotypePhenotype> m_diplotypes;

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

  public void addHaplotypeFunction(String haplotype, String func) {
    if (m_haplotypes != null) {
      m_haplotypes.put(haplotype, func);
    }
  }

  @Nullable
  public String lookupHaplotype(@Nullable String hap) {
    if (hap == null) {
      return null;
    }
    return m_haplotypes.get(hap);
  }

  /**
   * List of all diplotype to phenotype mappings for this gene
   */
  public List<DiplotypePhenotype> getDiplotypes() {
    return m_diplotypes;
  }

  public void setDiplotypes(List<DiplotypePhenotype> diplotypes) {
    m_diplotypes = diplotypes;
  }

  private String[] makePhenoPair(String diplotype) {
    String[] haps = diplotype.split("/");
    if (haps.length != 2) {
      throw new UnexpectedStateException("Diplotype doesn't have two alleles");
    }

    String hap1 = haps[0];
    String hap2 = haps[1];
    String func1 = getHaplotypes().get(hap1);
    String func2 = getHaplotypes().get(hap2);

    if (func1 == null) {
      sf_logger.warn("No function phenotype for " + getGene() + " " + hap1);
    } 
    if (func2 == null) {
      sf_logger.warn("No function phenotype for " + getGene() + " " + hap2);
    } 
    if (func1 != null && func2 != null) {
      return new String[]{func1,func2};
    } else {
      return null;
    }
  }

  public String makePhenotype(String diplotype) {

    String[] phenoPair = makePhenoPair(diplotype);
    if (phenoPair == null) {
      return "N/A";
    }

    Set<String> phenos = getDiplotypes().stream()
        .filter(d -> ((phenoPair[0].equals(d.getDiplotype().get(0)) && phenoPair[1].equals(d.getDiplotype().get(1))) || (phenoPair[0].equals(d.getDiplotype().get(1)) && phenoPair[1].equals(d.getDiplotype().get(0)))))
        .map(DiplotypePhenotype::getPhenotype)
        .collect(Collectors.toSet());
    if (phenos.size()>1) {
      throw new IllegalStateException("More than one phenotype match made");
    } else if (phenos.size() == 0) {
      return "N/A";
    } else {
      return phenos.iterator().next();
    }
  }

  public String toString() {
    return m_gene;
  }
}
