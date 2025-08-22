package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.common.util.ComparatorUtils;
import org.pharmgkb.pharmcat.reporter.TextConstants;


public class Genotype implements Comparable<Genotype> {
  @Expose
  @SerializedName("diplotypes")
  private final SortedSet<Diplotype> m_diplotypes = new TreeSet<>();

  private transient List<Map<String, Object>> m_lookupKeys;

  private transient SortedMap<String, Object> m_diplotypeKey = new TreeMap<>();


  /**
   * Private constructor for GSON.
   */
  @SuppressWarnings("unused")
  private Genotype() {
  }

  /**
   * Private constructor for tests.
   */
  private Genotype(Collection<Diplotype> diplotypes) {
    diplotypes.forEach(this::addDiplotype);
  }

  private Genotype(Diplotype diplotype) {
    addDiplotype(diplotype);
  }


  public SortedSet<Diplotype> getDiplotypes() {
    return m_diplotypes;
  }

  private void addDiplotype(Diplotype diplotype) {
    if (getDiplotypes().stream().anyMatch(d -> d.getGene().equals(diplotype.getGene()))) {
      throw new RuntimeException("Genotype already has diplotype for " + diplotype.getGene());
    }
    m_diplotypes.add(diplotype);

    m_diplotypeKey.put(diplotype.getGene(), diplotype.getDiplotypeKey());

    if (m_lookupKeys == null) {
      m_lookupKeys = diplotype.getLookupKeys().stream()
          .map((k) -> {
            Map<String, Object> lookupMap = new HashMap<>();
            lookupMap.put(diplotype.getGene(), k);
            return lookupMap;
          })
          .toList();
    } else {
      List<Map<String,Object>> newKeys = new ArrayList<>();
      for (String lookupKey : diplotype.getLookupKeys()) {
        for (Map<String, Object> originalKey : m_lookupKeys) {
          Map<String,Object> lookupMap = new HashMap<>(originalKey);
          lookupMap.put(diplotype.getGene(), lookupKey);
          newKeys.add(lookupMap);
        }
      }
      m_lookupKeys = newKeys;
    }
  }

  public boolean isInferred() {
    return m_diplotypes.stream().anyMatch(Diplotype::isInferred);
  }

  /**
   * Does any diplotype in this genotype rely on activity score for phenotype assignment?
   * @return true if any diplotype in this genotype uses AS
   */
  public boolean usesActivityScore() {
    return m_diplotypes.stream().anyMatch(Diplotype::isActivityScoreType);
  }

  public String toString() {
    if (m_diplotypes.isEmpty()) {
      return TextConstants.UNKNOWN_GENOTYPE;
    } else {
      return m_diplotypes.stream()
          .map(Diplotype::toString)
          .collect(Collectors.joining("; "));
    }
  }

  public String getPhenotypes() {
    if (m_diplotypes.isEmpty()) {
      return Haplotype.UNKNOWN;
    } else {
      return m_diplotypes.stream()
          .map(Diplotype::toString)
          .collect(Collectors.joining("; "));
    }
  }

  public List<Map<String, Object>> getLookupKeys() {
    return m_lookupKeys;
  }

  public Map<String, Object> getDiplotypeKey() {
    return m_diplotypeKey;
  }


  /**
   * Make a List of possible {@link Genotype} objects from different combinations of reporter diplotypes in the given
   * {@link GeneReport} objects.
   *
   * @param geneReports the {@link GeneReport} objects containing diplotypes to include in the possible genotypes
   * @return a List of all possible genotypes
   */
  public static List<Genotype> makeGenotypes(Collection<GeneReport> geneReports) {
    List<Genotype> possibleGenotypes = new ArrayList<>();

    for (GeneReport geneReport : geneReports) {
      if (possibleGenotypes.isEmpty()) {
        for (Diplotype diplotype : geneReport.getRecommendationDiplotypes()) {
          possibleGenotypes.add(new Genotype(diplotype));
        }
      } else {
        List<Genotype> oldGenotypes = possibleGenotypes;
        possibleGenotypes = new ArrayList<>();
        for (Diplotype diplotype : geneReport.getRecommendationDiplotypes()) {
          for (Genotype oldGenotype : oldGenotypes) {
            Genotype newGenotype = new Genotype(oldGenotype.getDiplotypes());
            newGenotype.addDiplotype(diplotype);
            possibleGenotypes.add(newGenotype);
          }
        }
      }
    }
    return possibleGenotypes;
  }


  public static Genotype forTest(List<Diplotype> diplotypes) {
    return new Genotype(diplotypes);
  }

  @Override
  public int compareTo(Genotype o) {
    if (o == this) {
      return 0;
    }
    return ComparatorUtils.compareCollection(m_diplotypes, o.getDiplotypes());
  }
}
