
package org.pharmgkb.pharmcat.haplotype.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.haplotype.MatchData;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;


/**
 * The {@link NamedAlleleMatcher} results for a single gene.
 *
 * @author Mark Woon
 */
public class GeneCall {
  @Expose
  @SerializedName("source")
  private final DataSource m_source;
  @Expose
  @SerializedName("version")
  private final String m_version;
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
  @SerializedName("haplotypeMatches")
  private final List<HaplotypeMatch> m_haplotypeMatches = new ArrayList<>();
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
  @Expose
  @SerializedName("warnings")
  private final List<MessageAnnotation> m_warnings;


  public GeneCall(DataSource source, String version, String chromosome, String gene,
      MatchData matchData, Set<String> uncallableHaplotypes, Set<String> ignoredHaplotypes,
      @Nullable List<MessageAnnotation> warnings) {
    m_source = source;
    m_version = version;
    m_chromosome = chromosome;
    m_gene = gene;
    m_matchData = matchData;
    m_uncallableHaplotypes = uncallableHaplotypes;
    m_ignoredHaplotypes = ignoredHaplotypes;
    m_variantsOfInterest = matchData.getExtraPositions();
    m_warnings = warnings;
  }


  /**
   * Gets the source of the definition file used to make this call.
   */
  public DataSource getSource() {
    return m_source;
  }

  /**
   * Gets the version of the definition file used to make this call.
   */
  public String getVersion() {
    return m_version;
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
    if (diplotype.getHaplotype2() != null) {
      m_haplotypes.add(diplotype.getHaplotype2());
    }
  }


  /**
   * Gets the haplotypes that are part of the called diplotypes.
   */
  public SortedSet<BaseMatch> getHaplotypes() {
    return m_haplotypes;
  }


  /**
   * Gets possible haplotype matches when no diplotypes are called.
   * <p>
   * This is currently only used for DPYD.
   */
  public List<HaplotypeMatch> getHaplotypeMatches() {
    return m_haplotypeMatches;
  }

  /**
   * Add haplotypes that were potential matches but could not be
   */
  public void addHaplotypeMatches(List<HaplotypeMatch> haplotypes) {
    m_haplotypeMatches.addAll(haplotypes);
    Collections.sort(m_haplotypeMatches);
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


  /**
   * Gets all warnings generated for this call.
   */
  public List<MessageAnnotation> getWarnings() {
    return m_warnings;
  }

  @Override
  public String toString() {
    return m_gene;
  }
}
