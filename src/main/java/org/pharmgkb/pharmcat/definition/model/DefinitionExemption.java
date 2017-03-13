package org.pharmgkb.pharmcat.definition.model;

import java.util.Collections;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;


/**
 * @author Mark Woon
 */
public class DefinitionExemption {
  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("ignoredAlleles")
  private SortedSet<String> m_ignoredAlleles;
  @Expose
  @SerializedName("ignoredAllelesLc")
  private SortedSet<String> m_ignoredAllelesLc;
  @Expose
  @SerializedName("allHits")
  private boolean m_allHits;
  @Expose
  @SerializedName("assumeReference")
  private boolean m_assumeReference;


  public DefinitionExemption(@Nonnull String gene, @Nullable SortedSet<String> ignoredAlleles, boolean allHits,
      boolean assumeReference) {
    m_gene = gene;
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
   * Get named alleles to ignore.
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
