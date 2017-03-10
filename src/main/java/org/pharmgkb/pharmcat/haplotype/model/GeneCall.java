
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.LinkedHashSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.haplotype.MatchData;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;


/**
 * The {@link NamedAlleleMatcher} results for a single gene.
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
  private LinkedHashSet<DiplotypeMatch> m_diplotypes = new LinkedHashSet<>();
  @Expose
  @SerializedName("haplotypes")
  private SortedSet<HaplotypeMatch> m_haplotypes = new TreeSet<>();
  @Expose
  @SerializedName("phased")
  private boolean m_isPhased = true;
  @Expose
  @SerializedName("variants")
  private SortedSet<Variant> m_variants = new TreeSet<>();
  @Expose
  @SerializedName("matchData")
  private MatchData m_matchData;
  @Expose
  @SerializedName("uncallableHaplotypes")
  private Set<String> m_uncallableHaplotypes;
  @Expose
  @SerializedName("ignoredHaplotypes")
  private Set<String> m_ignoredHaplotypes;


  public GeneCall(@Nonnull String alleleDefinitionVersion, @Nonnull String chromosome, @Nonnull String gene,
      @Nonnull MatchData matchData, @Nonnull Set<String> uncallableHaplotypes, @Nonnull Set<String> ignoredHaplotypes) {

    m_alleleDefinitionVersion = alleleDefinitionVersion;
    m_chromosome = chromosome;
    m_gene = gene;
    m_matchData = matchData;
    m_uncallableHaplotypes = uncallableHaplotypes;
    m_ignoredHaplotypes = ignoredHaplotypes;
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

  public @Nonnull Set<String> getIgnoredHaplotypes() {
    return m_ignoredHaplotypes;
  }


  public LinkedHashSet<DiplotypeMatch> getDiplotypes() {
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

  public void add(@Nonnull Variant pos) {
    Preconditions.checkNotNull(pos);
    m_variants.add(pos);
    if (!pos.isPhased()) {
      m_isPhased = false;
    }
  }


  /**
   * Gets if call was entirely based on phased data.
   * This will return false if any variants were unphased.
   */
  public boolean isPhased() {
    return m_isPhased;
  }


  @Override
  public String toString() {
    return m_gene;
  }
}
