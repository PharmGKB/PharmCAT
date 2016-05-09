
package org.pharmgkb.pharmcat.haplotype.model.json;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


public class GeneCall {
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
  private Set<DiplotypeMatch> diplotypes = new HashSet<>();
  @SerializedName("haplotypes")
  @Expose
  private SortedSet<HaplotypeMatch> m_haplotypes = new TreeSet<>();
  @SerializedName("variants")
  @Expose
  private SortedSet<Variant> m_variants = new TreeSet<>();
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


  public Set<DiplotypeMatch> getDiplotypes() {
    return diplotypes;
  }

  public void setDiplotypes(Set<DiplotypeMatch> diplotypes) {
    this.diplotypes = diplotypes;
  }

  public void addDiplotype(DiplotypeMatch diplotype) {
    diplotypes.add(diplotype);
    m_haplotypes.add(diplotype.getHaplotype1());
    m_haplotypes.add(diplotype.getHaplotype2());
  }


  public Set<HaplotypeMatch> getHaplotypes() {
    return m_haplotypes;
  }


  public SortedSet<Variant> getVariants() {
    return m_variants;
  }

  public void setVariants(SortedSet<Variant> variants) {
    m_variants = variants;
  }

  public void add(Variant pos) {
    m_variants.add(pos);
  }


  public List<HaplotypesNotCalled> getHaplotypesNotCalled() {
    return haplotypesNotCalled;
  }

  public void setHaplotypesNotCalled(List<HaplotypesNotCalled> haplotypesNotCalled) {
    this.haplotypesNotCalled = haplotypesNotCalled;
  }
}
