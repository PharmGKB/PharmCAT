package org.pharmgkb.pharmcat.phenotype.model;

import java.util.*;
import java.util.stream.Collectors;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.Constants;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * Model object that contains all the phenotype mapping data needed for a gene.
 *
 * @author Ryan Whaley
 */
public class GenePhenotype {
  public static final String UNASSIGNED_FUNCTION = "Unassigned function";

  @Expose
  @SerializedName("gene")
  private String m_gene;
  @Expose
  @SerializedName("haplotypes")
  private Map<String, String> m_haplotypes;
  @Expose
  @SerializedName("activityValues")
  private Map<String, String> m_activityValues = new HashMap<>();
  @Expose
  @SerializedName("diplotypes")
  private SortedSet<DiplotypeRecord> m_diplotypes;
  @Expose
  @SerializedName("namedAlleles")
  private List<HaplotypeRecord> m_namedAlleles;
  @Expose
  @SerializedName(value = "version")
  private String m_version;
  @Expose(serialize = false)
  @SerializedName("diplotypeFunctions")
  private List<DiplotypeFunction> m_diplotypeFunctions;

  // only used by Subsetter
  private final transient Set<String> m_modifiedActivity = new HashSet<>();
  private final transient Set<String> m_modifiedFunction = new HashSet<>();


  /**
   * The HGNC gene symbol
   */
  public String getGene() {
    return m_gene;
  }


  /**
   * Map of haplotype name (e.g. *1) to a function (e.g. Increased function)
   */
  public Map<String, String> getHaplotypes() {
    return m_haplotypes;
  }

  /**
   * Map of a haplotype name (e.g. *1) to an activity value (e.g. 1.0) if applicable. Not all genes use activity values
   */
  public Map<String,String> getActivityValues() {
    return m_activityValues;
  }

  private boolean isActivityGene() {
    return !m_activityValues.isEmpty();
  }

  /**
   * Makes a map combining the allele function and activity value (if applicable) and mapping it by allele name (e.g. *1)
   * @return a map of allele name to formatted function + score
   */
  public Map<String,String> makeFormattedFunctionScoreMap() {
    Map<String,String> newMap = new HashMap<>();
    for (HaplotypeRecord haplotypeRecord : getNamedAlleles()) {
      newMap.put(haplotypeRecord.getName(), haplotypeRecord.toFormattedFunction());
    }
    return newMap;
  }

  public void assignActivity(@Nullable Haplotype haplotype) {
    if (haplotype == null || haplotype.isUnknown()) {
      return;
    }
    haplotype.setActivityValue(m_activityValues.getOrDefault(haplotype.getName(), TextConstants.NA));
  }

  public boolean isMatchedByActivityScore() {
    return !m_activityValues.isEmpty() && m_activityValues.values().stream()
        .filter(StringUtils::isNotBlank)
        .anyMatch(v -> !v.equalsIgnoreCase(TextConstants.NA));
  }


  public String getHaplotypeFunction(String haplotype) {
    if (StringUtils.isBlank(haplotype) || m_namedAlleles == null) {
      return UNASSIGNED_FUNCTION;
    }

    return m_namedAlleles.stream()
        .filter(n -> n.getName().equals(haplotype))
        .findFirst()
        .map(HaplotypeRecord::getFunctionValue)
        .orElse(UNASSIGNED_FUNCTION);
  }

  public @Nullable String getHaplotypeActivity(String haplotype) {
    return m_activityValues.get(haplotype);
  }

  public @Nullable Float getHaplotypeActivityScore(String haplotype) {
    String value = m_activityValues.get(haplotype);
    if (value == null || value.equalsIgnoreCase(TextConstants.NA)) {
      return null;
    }
    return Float.valueOf(value);
  }


  /**
   * List of all diplotype to phenotype mappings for this gene
   */
  public SortedSet<DiplotypeRecord> getDiplotypes() {
    return m_diplotypes;
  }


  public List<HaplotypeRecord> getNamedAlleles() {
    return m_namedAlleles;
  }

  public Optional<DiplotypeRecord> findDiplotype(Map<String,Integer> diplotypeKey) {
    List<DiplotypeRecord> diplotypes = m_diplotypes.stream()
        .filter(d -> d.matchesKey(diplotypeKey))
        .toList();
    if (diplotypes.size() == 1) {
      return Optional.of(diplotypes.get(0));
    } else if (diplotypes.isEmpty()) {
      return Optional.empty();
    }
    // should never happen, DataManager should have caught this
    throw new IllegalStateException(diplotypes.size() + " diplotypes found for " + m_gene + " for " +
        diplotypeKey.keySet().stream()
            .map(k -> k + " (" + diplotypeKey.get(k) + ")")
            .collect(Collectors.joining(", ")));
  }


  /**
   * Gets the version of the {@link DataSource} the phenotype mapping is from.
   */
  public String getVersion() {
    return m_version;
  }


  /**
   * This is only available during data ingestion.
   */
  public List<DiplotypeFunction> getDiplotypeFunctions() {
    return this.m_diplotypeFunctions;
  }


  @Override
  public String toString() {
    return m_gene;
  }


  public void addHaplotypeRecord(String name, @Nullable String activityValue, @Nullable String functionValue,
      @Nullable String lookupKey) {
    if (isActivityGene()) {
      if (activityValue == null) {
        throw new IllegalStateException("Cannot add activity gene haplotype without activity value");
      }
      if (lookupKey == null) {
        lookupKey = activityValue;
      }
    } else {
      if (functionValue == null) {
        throw new IllegalStateException("Cannot add function gene haplotype without function value");
      }
      if (lookupKey == null) {
        lookupKey = functionValue;
      }
    }
    HaplotypeRecord hr = new HaplotypeRecord(name, activityValue, functionValue, lookupKey);
    m_namedAlleles.add(hr);
    m_haplotypes.put(name, hr.getLookupKey());
    if (hr.getActivityValue() != null) {
      m_activityValues.put(name, hr.getActivityValue());
    }
  }

  public boolean update(String allele, @Nullable String activityScore, @Nullable String function, DataSource src) {
    HaplotypeRecord hr = m_namedAlleles.stream()
        .filter(na -> na.getName().equals(allele))
        .findFirst()
        .orElse(null);
    if (hr == null) {
      // not all haplotypes will have a phenotype associated with it
      System.out.println("Adding phenotype for " + m_gene + " " + allele);
      addHaplotypeRecord(allele, activityScore, function, null);
      return true;
    }
    boolean gotChange = false;
    boolean modified = false;
    if (activityScore != null) {
      if (hr.getActivityValue() == null) {
        System.out.println("New " + src + " activity score for " + m_gene + " " + allele);
        modified = true;
      } else if (!hr.getActivityValue().equals(activityScore)){
        if (Double.parseDouble(hr.getActivityValue()) != Double.parseDouble(activityScore)) {
          System.out.println("Overwriting " + src + " activity score for " + m_gene + " " + allele +
              " (" + hr.getActivityValue() + " to " + activityScore + ")");
          modified = true;
        }
      }
    } else if (hr.getActivityValue() != null){
      System.out.println("Nulling out " + src + " activity score for " + m_gene + " " + allele);
      modified = true;
    }
    if (modified) {
      m_modifiedActivity.add(allele);
      hr.setActivityValue(activityScore);
      m_activityValues.put(allele, activityScore);
      gotChange = true;
    }

    modified = false;
    if (function != null) {
      if (hr.getFunctionValue() == null) {
        System.out.println("New " + src + " function for " + m_gene + " " + allele);
        modified = true;
      } else if (!hr.getFunctionValue().equals(function)) {
        System.out.println("Overwriting " + src + " function for " + m_gene + " " + allele + " (" +
            hr.getFunctionValue() + " to " + function + ")");
        modified = true;
      }
    } else if (hr.getFunctionValue() != null){
      System.out.println("Nulling out " + src + " function for " + m_gene + " " + allele);
      modified = true;
    }
    if (modified) {
      m_modifiedFunction.add(allele);
      hr.setFunctionValue(function);
      m_haplotypes.put(allele, function);
      gotChange = true;
    }
    return gotChange;
  }

  public boolean isModified(String allele) {
    return m_modifiedActivity.contains(allele) || m_modifiedFunction.contains(allele);
  }


  /**
   * Generates diplotypes based on {@link HaplotypeRecord}s.
   * <p>
   * NOT PART OF PUBLIC API.  Only used during data ingestion.
   */
  public void generateDiplotypes(DataSource source) {
    m_diplotypes = new TreeSet<>(makeDiplotypes(this, source));
  }

  private static Set<DiplotypeRecord> makeDiplotypes(GenePhenotype gp, DataSource source) {
    Set<DiplotypeRecord> results = new HashSet<>();
    Multimap<String,String> functionToAlleleMap = HashMultimap.create();
    for (HaplotypeRecord haplotype : gp.getNamedAlleles()) {
      functionToAlleleMap.put(haplotype.getLookupKey(), haplotype.getName());
    }

    for (DiplotypeFunction diplotypeFunction : gp.getDiplotypeFunctions()) {
      results.addAll(makeDiplotypes(diplotypeFunction, functionToAlleleMap, Constants.isActivityScoreGene(gp.getGene(), source)));
    }
    return results;
  }

  private static Set<DiplotypeRecord> makeDiplotypes(
      DiplotypeFunction diplotypeFunction,
      Multimap<String, String> functionToAllelesMap,
      boolean isActivityGene
  ) {
    Set<DiplotypeRecord> objects = new HashSet<>();
    Stack<String> fns = new Stack<>();
    String lookupKey = isActivityGene ? diplotypeFunction.activityScore : diplotypeFunction.phenotype;
    String activityScoreNormalized = !isActivityGene && TextConstants.isUnspecified(diplotypeFunction.activityScore) ? null : diplotypeFunction.activityScore;

    for (String key : diplotypeFunction.getLookupKey().keySet()) {
      for (int i = 0; i < diplotypeFunction.getLookupKey().get(key); i++) {
        fns.push(key);
      }
    }

    if (fns.size() == 2) {
      String fnKey1 = fns.pop();
      String fnKey2 = fns.pop();
      for (String allele1 : functionToAllelesMap.get(fnKey1)) {
        for (String allele2 : functionToAllelesMap.get(fnKey2)) {
          objects.add(makeDiplotype(diplotypeFunction.getPhenotype(), diplotypeFunction.getLookupKey(), allele1, allele2, activityScoreNormalized, lookupKey));
        }
      }
    } else {
      for (String allele : functionToAllelesMap.get(fns.pop())) {
        objects.add(makeDiplotype(diplotypeFunction.getPhenotype(), diplotypeFunction.getLookupKey(), allele, null, activityScoreNormalized, lookupKey));
      }
    }

    return objects;
  }

  private static DiplotypeRecord makeDiplotype(
      String phenotype, Map<String,Integer> phenotypeKey,
      String allele1, @Nullable String allele2, @Nullable String activityScore, String lookupKey
  ) {
    String diplotypeName;
    if (allele2 == null) {
      diplotypeName = allele1;

    } else {
      diplotypeName = Arrays.stream(new String[]{ allele1, allele2 })
          .sorted(HaplotypeNameComparator.getComparator())
          .collect(Collectors.joining("/"));
    }

    return new DiplotypeRecord(phenotype, diplotypeName, null, phenotype, makeDiplotypeKey(allele1, allele2),
        activityScore, lookupKey);
  }

  private static SortedMap<String, Integer> makeDiplotypeKey(String allele1, @Nullable String allele2) {
    SortedMap<String, Integer> diplotypeKey = new TreeMap<>(HaplotypeNameComparator.getComparator());
    if (Objects.equals(allele1, allele2)) {
      diplotypeKey.put(allele1, 2);
    } else {
      diplotypeKey.put(allele1, 1);
      if (allele2 != null) {
        diplotypeKey.put(allele2, 1);
      }
    }
    return diplotypeKey;
  }
}
