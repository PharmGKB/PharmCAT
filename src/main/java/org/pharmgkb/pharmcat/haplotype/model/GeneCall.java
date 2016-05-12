
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.HashSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.haplotype.MatchData;


/**
 * The Haplotyper results for a single gene.
 *
 * @author Mark Woon
 */
public class GeneCall {
  @Expose
  @SerializedName("alleleDefinitionVersion")
  private String m_alleleDefinitionVersion;
  @Expose
  @SerializedName("chromosome")
  private String m_chromosome;
  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("diplotypes")
  private Set<DiplotypeMatch> m_diplotypes = new HashSet<>();
  @Expose
  @SerializedName("haplotypes")
  private SortedSet<HaplotypeMatch> m_haplotypes = new TreeSet<>();
  @Expose
  @SerializedName("variants")
  private SortedSet<Variant> m_variants = new TreeSet<>();
  private MatchData m_matchData;
  private Set<String> m_uncallableHaplotypes;


  public GeneCall(@Nonnull String alleleDefinitionVersion, @Nonnull String chromosome, @Nonnull String gene,
      @Nonnull MatchData matchData, @Nonnull Set<String> uncallableHaplotypes) {

    m_alleleDefinitionVersion = alleleDefinitionVersion;
    m_chromosome = chromosome;
    m_gene = gene;
    m_matchData = matchData;
    m_uncallableHaplotypes = uncallableHaplotypes;
  }


  /**
   * Gets the version of the definition file used to make this call.
   */
  public @Nonnull String getAlleleDefinitionVersion() {
    return m_alleleDefinitionVersion;
  }

  public @Nonnull String getChromosome() {
    return m_chromosome;
  }

  public @Nonnull String getGene() {
    return m_gene;
  }

  public @Nonnull MatchData getMatchData() {
    return m_matchData;
  }

  public @Nonnull Set<String> getUncallableHaplotypes() {
    return m_uncallableHaplotypes;
  }


  public Set<DiplotypeMatch> getDiplotypes() {
    return m_diplotypes;
  }

  public void addDiplotype(DiplotypeMatch diplotype) {
    m_diplotypes.add(diplotype);
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
}
