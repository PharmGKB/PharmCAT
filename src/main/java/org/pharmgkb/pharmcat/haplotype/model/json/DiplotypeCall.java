
package org.pharmgkb.pharmcat.haplotype.model.json;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


public class DiplotypeCall {
  @SerializedName("gene")
  @Expose
  private String gene;
  @SerializedName("geneVersion")
  @Expose
  private String geneVersion;
  @SerializedName("chromosome")
  @Expose
  private String chromosome;
  @SerializedName("diplotypes")
  @Expose
  private Set<String> diplotypes = new HashSet<>();
  @SerializedName("haplotypes")
  @Expose
  private Set<HaplotypeCall> m_haplotypes = new HashSet<>();
  @SerializedName("variants")
  @Expose
  private List<Variant> m_variants = new ArrayList<>();
  @SerializedName("haplotypesNotCalled")
  @Expose
  private List<HaplotypesNotCalled> haplotypesNotCalled = new ArrayList<>();



  public String getGene() {
    return gene;
  }

  public void setGene(String gene) {
    this.gene = gene;
  }


  /**
   * Gets the version of the definition file used to make this call.
   */
  public String getGeneVersion() {
    return geneVersion;
  }

  public void setGeneVersion(String geneVersion) {
    this.geneVersion = geneVersion;
  }


  public String getChromosome() {
    return chromosome;
  }

  public void setChromosome(String chromosome) {
    this.chromosome = chromosome;
  }


  public Set<String> getDiplotypes() {
    return diplotypes;
  }

  public void setDiplotypes(Set<String> diplotypes) {
    this.diplotypes = diplotypes;
  }

  public void addDiplotype(String diplotype) {
    diplotypes.add(diplotype);
  }


  public Set<HaplotypeCall> getHaplotypes() {
    return m_haplotypes;
  }

  public void setHaplotypes(Set<HaplotypeCall> haplotypes) {
    m_haplotypes = haplotypes;
  }


  public List<Variant> getVariants() {
    return m_variants;
  }

  public void setVariants(List<Variant> variants) {
    this.m_variants = variants;
  }

  public void add(Variant pos) {
    this.m_variants.add(pos);
  }


  public List<HaplotypesNotCalled> getHaplotypesNotCalled() {
    return haplotypesNotCalled;
  }

  public void setHaplotypesNotCalled(List<HaplotypesNotCalled> haplotypesNotCalled) {
    this.haplotypesNotCalled = haplotypesNotCalled;
  }
}
