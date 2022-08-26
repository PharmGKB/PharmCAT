
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.Collection;
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
  @SerializedName("cpicVersion")
  private final String m_cpicVersion;
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
  private final SortedSet<BaseMatch> m_haplotypes = new TreeSet<>();
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


  public GeneCall(String alleleDefinitionVersion, String cpicVersion, String chromosome, String gene,
      MatchData matchData, Set<String> uncallableHaplotypes, Set<String> ignoredHaplotypes) {

    m_alleleDefinitionVersion = alleleDefinitionVersion;
    m_cpicVersion = cpicVersion;
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

  public String getCpicVersion() {
    return m_cpicVersion;
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

  public void addAllHaplotypes(Collection<HaplotypeMatch> haplotypes) {
    m_haplotypes.addAll(haplotypes);
  }


  public Set<BaseMatch> getHaplotypes() {
    return m_haplotypes;
  }


  public SortedSet<Variant> getVariants() {
    return m_variants;
  }

  /**
   * Adds a {@link Variant} to this {@link GeneCall}.
   * <p>
   * Call {@link #finalizeVariants()} when done adding.
   */
  public void add(Variant pos) {
    Preconditions.checkNotNull(pos);
    m_variants.add(pos);
    if (!pos.isPhased()) {
      m_isPhased = false;
    }
  }

  /**
   * Call this when done adding {@link Variant}s to make sure phasing is set correctly.
   */
  public void finalizeVariants() {
    if (m_variants.isEmpty()) {
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

  /**
   * Gets if call was based on "effectively phased" data (i.e. actually phased or unphased but homozygous at all
   * positions).
   */
  public boolean isEffectivelyPhased() {
    return m_matchData.isEffectivelyPhased();
  }


  @Override
  public String toString() {
    return m_gene;
  }
}
