package org.pharmgkb.pharmcat.phenotype.model;

import java.util.Map;
import com.google.common.base.Objects;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.ComparisonChain;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * Model that represents the diplotype data in the PharmGKB database.
 *
 * @author Ryan Whaley
 */
public class DiplotypeRecord implements Comparable<DiplotypeRecord> {

  @SerializedName("generesult")
  @Expose
  private String m_geneResult;
  @SerializedName("diplotype")
  @Expose
  private String m_diplotype;
  @SerializedName("description")
  @Expose
  private @Nullable String m_description;
  @SerializedName("lookupkey")
  @Expose
  private String m_lookupKey;
  @SerializedName("diplotypekey")
  @Expose
  private Map<String,Integer> m_diplotypeKey;
  @SerializedName("activityScore")
  @Expose
  private @Nullable String m_activityScore;
  @SerializedName("phenotype")
  @Expose
  private String m_phenotype;


  public DiplotypeRecord(String geneResult, String diplotype, @Nullable String description, String phenotype,
      Map<String,Integer> diplotypeKey, @Nullable String activityScore, String lookupKey) {
    this.m_phenotype = phenotype;
    this.m_geneResult = geneResult;

    this.m_diplotype = diplotype;
    this.m_description = description;
    this.m_diplotypeKey = diplotypeKey;
    this.m_activityScore = activityScore;
    this.m_lookupKey = lookupKey;
  }


  public String getGeneResult() {
    return m_geneResult;
  }

  /**
   * The bare diplotype name.
   *
   * @return diplotype in the form "*1/*3"
   */
  public String getDiplotype() {
    return m_diplotype;
  }

  public String toString() {
    return m_diplotype;
  }

  public @Nullable String getDescription() {
    return m_description;
  }

  public String getLookupKey() {
    return m_lookupKey;
  }

  public Map<String, Integer> getDiplotypeKey() {
    return m_diplotypeKey;
  }

  public @Nullable String getActivityScore() {
    return this.m_activityScore;
  }

  public String getPhenotype() {
    return this.m_phenotype;
  }


  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (obj == this) {
      return true;
    }
    if (obj.getClass() != getClass()) {
      return false;
    }
    DiplotypeRecord other = (DiplotypeRecord)obj;
    return Objects.equal(m_geneResult, other.getGeneResult()) &&
        Objects.equal(m_diplotype, other.getDiplotype()) &&
        Objects.equal(m_description, other.getDescription()) &&
        Objects.equal(m_lookupKey, other.getLookupKey()) &&
        Objects.equal(m_diplotypeKey, other.getDiplotypeKey());
  }

  public boolean matchesKey(Map<String,Integer> otherKey) {
    if (otherKey == null || otherKey.isEmpty()) {
      return false;
    }
    if (m_diplotypeKey == null || m_diplotypeKey.isEmpty()) {
      return false;
    }
    return otherKey.keySet().size() == m_diplotypeKey.keySet().size() && otherKey.entrySet().containsAll(m_diplotypeKey.entrySet());
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(m_geneResult, m_diplotype, m_description, m_lookupKey, m_diplotypeKey);
  }

  @Override
  public int compareTo(DiplotypeRecord o) {
    if (this == o) {
      return 0;
    }

    String[] dips1 = m_diplotype.split("/");
    String[] dips2 = o.getDiplotype().split("/");

    int rez = dips1.length - dips2.length;
    if (rez != 0) {
      return rez;
    }

    rez = HaplotypeNameComparator.getComparator().compare(dips1[0], dips2[0]);
    if (rez != 0) {
      return rez;
    }
    if (dips1.length > 1) {
      if (dips2.length > 1) {
        rez = HaplotypeNameComparator.getComparator().compare(dips1[1], dips2[1]);
        if (rez != 0) {
          return rez;
        }
      } else {
        return 1;
      }
    }
    return new ComparisonChain()
        .compare(m_activityScore, o.getActivityScore())
        .compare(m_geneResult, o.getGeneResult())
        .compare(m_phenotype, o.getPhenotype())
        .result();
  }
}
