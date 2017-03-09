package org.pharmgkb.pharmcat.definition.model;

import java.util.Collections;
import java.util.SortedSet;
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
  @SerializedName("allHits")
  private boolean m_allHits;


  public DefinitionExemption(@Nonnull String gene, @Nullable SortedSet<String> ignoredAlleles, boolean allHits) {
    m_gene = gene;
    if (ignoredAlleles == null) {
      m_ignoredAlleles = Collections.emptySortedSet();
    } else {
      m_ignoredAlleles = ignoredAlleles;
    }
    m_allHits = allHits;
  }


  public String getGene() {
    return m_gene;
  }

  public SortedSet<String> getIgnoredAlleles() {
    return m_ignoredAlleles;
  }

  public boolean isAllHits() {
    return m_allHits;
  }
}
