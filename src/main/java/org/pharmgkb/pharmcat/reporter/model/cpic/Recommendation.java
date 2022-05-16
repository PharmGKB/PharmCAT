package org.pharmgkb.pharmcat.reporter.model.cpic;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Predicate;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;


/**
 * A Recommendation object sourced from the CPIC DB. This will reflect the properties of the CPIC Recommendation data
 * model.
 */
public class Recommendation {
  @Expose
  @SerializedName("implications")
  private Map<String,String> m_implications;
  @Expose
  @SerializedName("drugrecommendation")
  private String m_drugRecommendation;
  @Expose
  @SerializedName("classification")
  private String m_classification;
  @Expose
  @SerializedName("phenotypes")
  private Map<String,String> m_phenotypes;
  @Expose
  @SerializedName("activityscore")
  private Map<String,String> m_activityScore;
  @Expose
  @SerializedName("allelestatus")
  private Map<String,String> m_alleleStatus;
  @Expose
  @SerializedName("lookupkey")
  private Map<String,String> m_lookupKey;
  @Expose
  @SerializedName("comments")
  private String m_comments;
  @Expose
  @SerializedName("population")
  private String m_population;
  @Expose
  @SerializedName("genotypes")
  private List<Map<String,String>> m_genotypes;
  @Expose
  @SerializedName("source")
  private DataSource m_source;

  private final Set<String> m_matchedDiplotypes = new TreeSet<>();
  private final List<Genotype> f_matchedGenotypes = new ArrayList<>();

  public Map<String, String> getImplications() {
    return m_implications;
  }

  public void setImplications(Map<String, String> implications) {
    m_implications = implications;
  }

  public String getDrugRecommendation() {
    return m_drugRecommendation;
  }

  public void setDrugRecommendation(String drugRecommendation) {
    m_drugRecommendation = drugRecommendation;
  }

  public String getClassification() {
    return m_classification;
  }

  public void setClassification(String classification) {
    m_classification = classification;
  }

  public Map<String, String> getPhenotypes() {
    return m_phenotypes;
  }

  public void setPhenotypes(Map<String, String> phenotypes) {
    m_phenotypes = phenotypes;
  }

  public Map<String, String> getActivityScore() {
    return m_activityScore;
  }

  public void setActivityScore(Map<String, String> activityScore) {
    m_activityScore = activityScore;
  }

  public Map<String, String> getAlleleStatus() {
    return m_alleleStatus;
  }

  public void setAlleleStatus(Map<String, String> alleleStatus) {
    m_alleleStatus = alleleStatus;
  }

  public Map<String, String> getLookupKey() {
    return m_lookupKey;
  }

  public void setLookupKey(Map<String, String> lookupKey) {
    m_lookupKey = lookupKey;
  }

  public Set<String> getMatchedDiplotypes() {
    return m_matchedDiplotypes;
  }

  public void addMatchedDiplotype(String dip) {
    m_matchedDiplotypes.add(dip);
  }

  public String getComments() {
    return m_comments;
  }

  public void setComments(String comments) {
    m_comments = comments;
  }

  public String getPopulation() {
    return m_population;
  }

  public void setPopulation(String population) {
    m_population = population;
  }

  public List<Map<String,String>> getGenotypes() {
    return m_genotypes;
  }

  protected void setGenotypes(List<Map<String,String>> genotypes) {
    m_genotypes = genotypes;
  }

  public DataSource getSource() {
    return m_source;
  }

  public void setSource(DataSource source) {
    m_source = source;
  }

  /**
   * Whether the given diplotypes match this Recommendation
   * @param diplotypes a collection of "GENE:*1/*4" diplotype strings, one per gene
   * @return true if a match
   */
  public boolean matchesGenotype(Collection<String> diplotypes) {
    if (diplotypes == null || diplotypes.size() == 0) {
      return false;
    }
    Map<String,String> parsedInput = new TreeMap<>();
    diplotypes.stream()
        .map(d -> d.split(":"))
        .forEach(tokens -> parsedInput.put(tokens[0], tokens[1]));

    return getGenotypes().stream()
        .anyMatch(g -> (parsedInput.keySet().containsAll(g.keySet())) && (
            g.keySet().stream().allMatch(key -> Objects.equals(g.get(key), parsedInput.get(key)))
            ));
  }

  public boolean matchesGenotype(Genotype genotype) {
    List<Map<String,String>> keys = genotype.toLookupKeys();
    return keys.stream()
        .anyMatch(matchesLookupKey);
  }

  protected final Predicate<Map<String,String>> matchesLookupKey = (lookupKey) -> getGenotypes() != null
      && getGenotypes().stream()
      .anyMatch(g -> (lookupKey.keySet().containsAll(g.keySet())) && (
          g.keySet().stream().allMatch(key -> Objects.equals(g.get(key), lookupKey.get(key)))
      ));

  public List<Genotype> getMatchedGenotypes() {
    return f_matchedGenotypes;
  }

  public void addMatchedGenotype(Genotype genotype) {
    f_matchedGenotypes.add(genotype);
  }
}
