package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.reflect.Type;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.Spliterator;
import java.util.TreeSet;
import java.util.function.Consumer;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.gson.reflect.TypeToken;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.DataSerializer;


/**
 * This class will read the drug data files from the filesystem and serve the data as {@link Drug} objects. This
 * includes data from all drug data sources.
 */
public class DrugCollection implements Iterable<Drug> {
  public static final String CPIC_FILE_NAME = "drugs.json";
  public static final String CPIC_URL = "https://files.cpicpgx.org/data/report/current/" + CPIC_FILE_NAME;
  public static final Type DRUG_LIST_TYPE = new TypeToken<ArrayList<Drug>>(){}.getType();
  public static final Path GUIDELINES_DIR =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/guidelines/cpic");

  private final Set<Drug> m_drugList = new TreeSet<>();
  private SortedSet<String> m_genes;


  /**
   * Default constructor.
   */
  public DrugCollection() throws IOException {
    this(GUIDELINES_DIR);
  }

  public DrugCollection(Path dir) throws IOException {
    List<Path> annotationFiles = new ArrayList<>();
    try (Stream<Path> stream = Files.list(dir)) { {
      stream.filter(f -> f.getFileName().toString().endsWith(".json"))
          .forEach(annotationFiles::add);
    }}
    if (annotationFiles.size() == 0) {
      throw new IOException("Cannot find annotations");
    }
    for (Path guidelineFile : annotationFiles) {
      try (BufferedReader br = Files.newBufferedReader(guidelineFile)) {
        m_drugList.add(DataSerializer.GSON.fromJson(br, Drug.class));
      }
    }
  }


  // A drug is ignored if it associated with ANY ignored gene.
  private static final Predicate<Drug> sf_filterIgnoredGenes = (drug) -> drug.getGenes().stream().noneMatch(GeneReport::isIgnored);

  public int size() {
    return m_drugList.size();
  }

  /**
   * Get all {@link Drug} objects
   * @return a List of {@link Drug} objects
   */
  public Set<Drug> list() {
    return m_drugList;
  }

  /**
   * Get only the "reportable" drugs that do not use "ignored" genes. A reportable drug cannot use _ANY_ ignored gene.
   * @return a List of {@link Drug} that can be reported on
   */
  public List<Drug> listReportable() {
    return m_drugList.stream().filter(sf_filterIgnoredGenes).collect(Collectors.toList());
  }

  /**
   * Find a {@link Drug} object that has the given ID or name
   * @return an optional {@link Drug} object
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
  public SortedSet<String> getAllReportableGenes() {
    if (m_genes == null) {
      m_genes = m_drugList.stream()
          .flatMap(d -> d.getGenes().stream())
          .filter(g -> !GeneReport.isIgnored(g))
          .collect(Collectors.toCollection(TreeSet::new));
    }
    return m_genes;
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

  public List<String> diff(DrugCollection otherCollection) {
    List<String> messages = new ArrayList<>();

    Set<String> otherDrugNames = otherCollection.list().stream()
        .map(Drug::getDrugName).collect(Collectors.toSet());

    Set<String> thisDrugNames = list().stream()
        .map(Drug::getDrugName).collect(Collectors.toSet());

    Set<String> removedDrugs = new HashSet<>(otherDrugNames);
    removedDrugs.removeAll(thisDrugNames);
    removedDrugs.stream().map(d -> String.format("Removed drug: %s", d)).forEach(messages::add);

    Set<String> addedDrugs = new HashSet<>(thisDrugNames);
    addedDrugs.removeAll(otherDrugNames);
    addedDrugs.stream().map(d -> String.format("Added drug: %s", d)).forEach(messages::add);

    return messages;
  }
}
