package org.pharmgkb.pharmcat.definition.model;

import java.util.Collections;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * @author Mark Woon
 */
public class DefinitionExemption {
  @Expose
  @SerializedName("gene")
  private final String m_gene;
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
  private final boolean m_allHits;
  @Expose
  @SerializedName("assumeReference")
  private final boolean m_assumeReference;


  public DefinitionExemption(String gene, @Nullable SortedSet<VariantLocus> extraPositions,
      @Nullable SortedSet<String> ignoredAlleles, boolean allHits, boolean assumeReference) {
    m_gene = gene;
    if (extraPositions == null) {
      m_extraPositions = Collections.emptySortedSet();
    } else {
      m_extraPositions = extraPositions;
    }
    if (ignoredAlleles == null) {
      m_ignoredAlleles = Collections.emptySortedSet();
      m_ignoredAllelesLc = m_ignoredAlleles;
    } else {
      m_ignoredAlleles = ignoredAlleles;
      m_ignoredAllelesLc = ignoredAlleles.stream().map(String::toLowerCase).collect(Collectors.toCollection(TreeSet::new));
    }
    m_allHits = allHits;
    m_assumeReference = assumeReference;
  }


  public String getGene() {
    return m_gene;
  }


  /**
   * Gets the extra positions to pull allele information for.
   */
  public SortedSet<VariantLocus> getExtraPositions() {
    return m_extraPositions;
  }

  /**
   * Gets the named alleles to ignore.
   */
  public SortedSet<String> getIgnoredAlleles() {
    return m_ignoredAlleles;
  }

  /**
   * Checks if should ignore the given named allele.
   */
  public boolean shouldIgnore(String allele) {
    return m_ignoredAllelesLc.contains(allele.toLowerCase());
  }

  /**
   * Gets if all possible diplotypes should be reported.
   * Default is false.
   */
  public boolean isAllHits() {
    return m_allHits;
  }


  /**
   * Gets if missing positions should be assumed to be reference.
   * Default is true.
   */
  public boolean isAssumeReference() {
    return m_assumeReference;
  }
}
