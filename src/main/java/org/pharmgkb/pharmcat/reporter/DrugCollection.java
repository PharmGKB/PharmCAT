package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Type;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.Spliterator;
import java.util.TreeSet;
import java.util.function.Consumer;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


/**
 * This class will read the drugs.json file and serve the data as {@link Drug} objects.
 */
public class DrugCollection implements Iterable<Drug> {
  private static final String CPIC_FILE_NAME = "drugs.json";
  private static final String DPWG_FILE_NAME = "drugsDpwg.json";
  private static final Type DRUG_LIST_TYPE = new TypeToken<ArrayList<Drug>>(){}.getType();
  private static final Gson GSON = new GsonBuilder()
      .serializeNulls()
      .excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting().create();

  private final Set<Drug> m_drugList = new TreeSet<>();

  /**
   * Default constructor. Will use the drugs.json file defined in the codebase
   * @throws IOException when the drugs.json file cannot be read
   */
  public DrugCollection() throws IOException {
    try (InputStream inputStream = getClass().getResourceAsStream(CPIC_FILE_NAME)) {
      if (inputStream == null) {
        throw new RuntimeException("Drug definition file not found");
      }
      addDrugsFromStream(inputStream, DataSource.CPIC);
    }

    try (InputStream inputStream = getClass().getResourceAsStream(DPWG_FILE_NAME)) {
      if (inputStream == null) {
        throw new RuntimeException("DPWG Drug definition file not found");
      }
      addDrugsFromStream(inputStream, DataSource.DPWG);
    }
  }

  /**
   * Alternate constructor. Will use the supplied stream to load drug data. Expected to deserialize into {@link Drug}
   * objects.
   * @param inputStream an InputStream of serialized Drug JSON data
   * @throws IOException when the input stream cannot be read
   */
  public DrugCollection(InputStream inputStream, DataSource source) throws IOException {
    addDrugsFromStream(inputStream, source);
  }

  private void addDrugsFromStream(InputStream inputStream, DataSource source) throws IOException {
    try (BufferedReader br = new BufferedReader(new InputStreamReader(inputStream))) {
      List<Drug> drugs = GSON.fromJson(br, DRUG_LIST_TYPE);
      drugs.forEach(d -> d.setSource(source));
      m_drugList.addAll(drugs);
    }
  }

  // A drug is ignored if it associated with ANY ignored gene.
  private static final Predicate<Drug> sf_filterIgnoredGenes = (drug) -> drug.getGenes().stream().noneMatch(GeneReport::isIgnored);
  public static final Predicate<Drug> ONLY_CPIC = (drug) -> drug.getSource() == DataSource.CPIC;
  public static final Predicate<Drug> ONLY_DPWG = (drug) -> drug.getSource() == DataSource.DPWG;

  public int size() {
    return m_drugList.size();
  }

  /**
   * Get all CPIC {@link Drug} objects
   * @return a List of CPIC {@link Drug} objects
   */
  public Set<Drug> list() {
    return m_drugList;
  }

  /**
   * Get only the "reportable" drugs that do note use "ignored" genes. A reportable drug cannot use _ANY_ ignored gene.
   * @return a List of {@link Drug} that can be reported on
   */
  public List<Drug> listReportable() {
    return m_drugList.stream().filter(sf_filterIgnoredGenes).collect(Collectors.toList());
  }

  /**
   * Find a CPIC {@link Drug} object that has the given ID or name
   * @return an optional CPIC {@link Drug} object
   */
  public Optional<Drug> find(String identifier) {
    if (StringUtils.isBlank(identifier)) return Optional.empty();
    return m_drugList.stream()
        .filter((d) -> d.getDrugId().equalsIgnoreCase(identifier) || d.getDrugName().equalsIgnoreCase(identifier))
        .findFirst();
  }

  /**
   * Gets a Set of all genes that are related to drugs for this collection and not ignored
   * @return a Set of all genes reportable by these drugs
   */
  public Set<String> getAllReportableGenes() {
    return m_drugList.stream()
        .flatMap(d -> d.getGenes().stream())
        .filter(g -> !GeneReport.isIgnored(g))
        .collect(Collectors.toSet());
  }

  @Override
  public Iterator<Drug> iterator() {
    return m_drugList.iterator();
  }

  @Override
  public void forEach(Consumer<? super Drug> action) {
    m_drugList.forEach(action);
  }

  @Override
  public Spliterator<Drug> spliterator() {
    return m_drugList.spliterator();
  }
}
