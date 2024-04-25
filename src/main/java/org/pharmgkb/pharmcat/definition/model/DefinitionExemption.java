package org.pharmgkb.pharmcat.definition.model;

import java.util.Collections;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * This class represents special exemptions that must be applied to the {@link DefinitionFile} for the given
 * {@code gene}.
 *
 * @author Mark Woon
 */
public class DefinitionExemption implements Comparable<DefinitionExemption> {
  @Expose
  @SerializedName("gene")
  private final String m_gene;
  @Expose
  @SerializedName("ignoredPositions")
  private final SortedSet<VariantLocus> m_ignoredPositions;
  @Expose
  @SerializedName("extraPositions")
  private final SortedSet<VariantLocus> m_extraPositions;
  @Expose
  @SerializedName("ignoredAlleles")
  private final SortedSet<String> m_ignoredAlleles;
  @Expose
  @SerializedName("ignoredAllelesLc")
  private final SortedSet<String> m_ignoredAllelesLc;
  @Expose
  @SerializedName("allHits")
  private final Boolean m_allHits;


  public DefinitionExemption(String gene, @Nullable SortedSet<VariantLocus> ignoredPositions,
      @Nullable SortedSet<VariantLocus> extraPositions, @Nullable SortedSet<String> ignoredAlleles, Boolean allHits) {
    m_gene = gene;

    m_ignoredPositions = Objects.requireNonNullElse(ignoredPositions, Collections.emptySortedSet());
    m_extraPositions = Objects.requireNonNullElse(extraPositions, Collections.emptySortedSet());
    if (ignoredAlleles == null) {
      m_ignoredAlleles = Collections.emptySortedSet();
      m_ignoredAllelesLc = m_ignoredAlleles;
    } else {
      m_ignoredAlleles = ignoredAlleles;
      m_ignoredAllelesLc = ignoredAlleles.stream()
          .map(String::toLowerCase)
          .collect(Collectors.toCollection(TreeSet::new));
    }
    m_allHits = allHits;
  }


  public String getGene() {
    return m_gene;
  }


  /**
   * Gets the positions from original definition that should be ignored.
   * These get removed by the {@link org.pharmgkb.pharmcat.util.DataManager} when definitions are first pulled back from
   * PharmGKB.
   */
  public SortedSet<VariantLocus> getIgnoredPositions() {
    return m_ignoredPositions;
  }

  /**
   * Checks if given position should be ignored.
   * <p>
   * <b>Currently only checks based on RSID!</b>
   */
  boolean shouldIgnorePosition(VariantLocus position) {
    return m_ignoredPositions.stream()
        .anyMatch(vl -> {
          if (vl.getRsid() != null) {
            return vl.getRsid().equals(position.getRsid());
          }
          return vl.getChromosomeHgvsName().equals(position.getChromosomeHgvsName());
        });
  }


  /**
   * Gets the extra positions to pull allele information for.
   */
  public SortedSet<VariantLocus> getExtraPositions() {
    return m_extraPositions;
  }


  /**
   * Gets the named alleles from the original definition that should be ignored.
   * These get removed by the {@link org.pharmgkb.pharmcat.util.DataManager} when definitions are first pulled back
   * from PharmGKB.
   */
  public SortedSet<String> getIgnoredAlleles() {
    return m_ignoredAlleles;
  }

  /**
   * Checks if the given named allele should be ignored.
   */
  boolean shouldIgnoreAllele(String allele) {
    return m_ignoredAllelesLc.contains(allele.toLowerCase());
  }

  /**
   * Gets if all possible diplotypes should be reported.
   */
  public @Nullable Boolean isAllHits() {
    return m_allHits;
  }


  @Override
  public int compareTo(DefinitionExemption o) {
    Preconditions.checkNotNull(m_gene);
    Preconditions.checkNotNull(o.getGene());
    return m_gene.compareTo(o.getGene());
  }
}
