
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.LinkedHashSet;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
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
  private final String m_alleleDefinitionVersion;
  @Expose
  @SerializedName("chromosome")
  private final String m_chromosome;
  @Expose
  @SerializedName("gene")
  private final String m_gene;
  @Expose
  @SerializedName("diplotypes")
  private final LinkedHashSet<DiplotypeMatch> m_diplotypes = new LinkedHashSet<>();
  @Expose
  @SerializedName("haplotypes")
  private final SortedSet<HaplotypeMatch> m_haplotypes = new TreeSet<>();
  @Expose
  @SerializedName("phased")
  private boolean m_isPhased = true;
  @Expose
  @SerializedName("variants")
  private SortedSet<Variant> m_variants = new TreeSet<>();
  @Expose
  @SerializedName("variantsOfInterest")
  private final SortedSet<Variant> m_variantsOfInterest;
  @Expose
  @SerializedName("matchData")
  private final MatchData m_matchData;
  @Expose
  @SerializedName("uncallableHaplotypes")
  private final Set<String> m_uncallableHaplotypes;
  @Expose
  @SerializedName("ignoredHaplotypes")
  private final Set<String> m_ignoredHaplotypes;


  public GeneCall(String alleleDefinitionVersion, String chromosome, String gene,
      MatchData matchData, Set<String> uncallableHaplotypes, Set<String> ignoredHaplotypes) {

    m_alleleDefinitionVersion = alleleDefinitionVersion;
    m_chromosome = chromosome;
    m_gene = gene;
    m_matchData = matchData;
    m_uncallableHaplotypes = uncallableHaplotypes;
    m_ignoredHaplotypes = ignoredHaplotypes;
    m_variantsOfInterest = matchData.getExtraPositions();
  }


  /**
   * Gets the version of the definition file used to make this call.
   */
  public String getAlleleDefinitionVersion() {
    return m_alleleDefinitionVersion;
  }

  public String getChromosome() {
    return m_chromosome;
  }

  public String getGene() {
    return m_gene;
  }

  public MatchData getMatchData() {
    return m_matchData;
  }

  public Set<String> getUncallableHaplotypes() {
    return m_uncallableHaplotypes;
  }

  public Set<String> getIgnoredHaplotypes() {
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

  public void add(Variant pos) {
    Preconditions.checkNotNull(pos);
    m_variants.add(pos);
    if (!pos.isPhased()) {
      m_isPhased = false;
    }
  }


  public SortedSet<Variant> getVariantsOfInterest() {
    return m_variantsOfInterest;
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
