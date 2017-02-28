package org.pharmgkb.pharmcat.definition.model;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.UnexpectedStateException;


/**
 * Model object that contains all the phenotype mapping data needed for a gene
 *
 * @author Ryan Whaley
 */
public class GenePhenotype {

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

    return new String[]{func1,func2};
  }

  public String makeFunction(String diplotype) {

    String[] phenoPair = makePhenoPair(diplotype);

    if (StringUtils.isNotBlank(phenoPair[0]) && StringUtils.isNotBlank(phenoPair[1])) {

      if (phenoPair[0].equals(phenoPair[1])) {
        return "Two " + phenoPair[0].toLowerCase() + " alleles";
      }
      else {
        return "One " + phenoPair[0].toLowerCase() + " allele and one " + phenoPair[1].toLowerCase() + " allele";
      }

    }
    return "N/A";
  }

  public String makePhenotype(String diplotype) {

    String[] phenoPair = makePhenoPair(diplotype);

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

}
