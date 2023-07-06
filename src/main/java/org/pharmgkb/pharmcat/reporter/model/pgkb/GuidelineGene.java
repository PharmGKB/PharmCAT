package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


/**
 * PharmGKB Guideline Annotation Gene Model. These genes are the ones that are involved in annotation group matching.
 */
@Deprecated
public class GuidelineGene {

  @SerializedName("alleles")
  @Expose
  private List<GuidelineAllele> m_alleles = new ArrayList<>();
  @SerializedName("gene")
  @Expose
  private Gene m_gene;

  public List<GuidelineAllele> getAlleles() {
    return m_alleles;
  }

  public void setAlleles(List<GuidelineAllele> alleles) {
    m_alleles = alleles;
  }

  public Gene getGene() {
    return m_gene;
  }

  public void setGene(Gene gene) {
    m_gene = gene;
  }

  private Optional<String> findFunctionForAllele(String alleleName) {
    return m_alleles.stream()
        .filter(a -> a.getLabel().equals(alleleName))
        .findFirst()
        .map(a -> a.getFunctionTerm().getTerm());
  }

  public Optional<String> findFunctionForAllele(Haplotype haplotype) {
    if (haplotype == null) {
      return Optional.empty();
    }
    return findFunctionForAllele(haplotype.getName());
  }
}
