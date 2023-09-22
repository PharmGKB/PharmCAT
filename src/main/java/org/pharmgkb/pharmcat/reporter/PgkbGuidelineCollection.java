package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.pgkb.AccessionObject;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;
import org.pharmgkb.pharmcat.util.DataSerializer;


public class PgkbGuidelineCollection {
  private static final Path GUIDELINES_DIR =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/guidelines");

  private final List<GuidelinePackage> f_guidelinePackages = new ArrayList<>();
  private final SortedSetMultimap<String,GuidelinePackage> f_guidelineMap = TreeMultimap.create(String::compareToIgnoreCase, Comparator.naturalOrder());
  private SortedSet<String> m_genes;


  public PgkbGuidelineCollection() throws IOException {
    this(GUIDELINES_DIR);
  }

  public PgkbGuidelineCollection(Path dir) throws IOException {
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
        GuidelinePackage guidelinePackage = DataSerializer.GSON.fromJson(br, GuidelinePackage.class);
        f_guidelinePackages.add(guidelinePackage);
        for (AccessionObject chemical : guidelinePackage.getGuideline().getRelatedChemicals()) {
          f_guidelineMap.put(chemical.getName(), guidelinePackage);
        }
      }
    }
  }

  public List<GuidelinePackage> getGuidelinePackages() {
    return f_guidelinePackages;
  }

  public List<GuidelinePackage> findGuidelinePackages(String chemicalName, DataSource source) {
    return f_guidelineMap.get(chemicalName).stream()
        .filter(p -> p.getGuideline().getSource().equalsIgnoreCase(source.getPharmgkbName()))
        .collect(Collectors.toList());
  }

  public SortedSetMultimap<String,GuidelinePackage> getGuidelineMap() {
    return f_guidelineMap;
  }

  public Set<GuidelinePackage> getGuidelinesFromSource(DataSource dataSource) {
    return f_guidelineMap.values().stream()
        .filter(g -> g.getGuideline().getSource().equalsIgnoreCase(dataSource.getPharmgkbName()))
        .collect(Collectors.toSet());
  }

  public Set<String> getChemicalsUsedInSource(DataSource source) {
    return f_guidelineMap.values().stream()
        .map(GuidelinePackage::getGuideline)
        .filter(g -> g.getSource().equalsIgnoreCase(source.getPharmgkbName()))
        .flatMap(g -> g.getRelatedChemicals().stream())
        .map(AccessionObject::getName)
        .collect(Collectors.toSet());
  }

  public SortedSet<String> getGenesWithRecommendations() {
    if (m_genes == null) {
      m_genes = f_guidelinePackages.stream()
          .flatMap(p -> p.getRecommendations().stream())
          .filter(Objects::nonNull)
          .flatMap(r -> r.getLookupKey().keySet().stream())
          .collect(Collectors.toCollection(TreeSet::new));
    }
    return m_genes;
  }

  public SortedSet<String> getGenesUsedInSource(DataSource source) {
    return f_guidelinePackages.stream()
        .filter(p -> p.getGuideline().getSource().equalsIgnoreCase(source.getPharmgkbName()))
        .flatMap(p -> p.getGenes().stream())
        .filter(Objects::nonNull)
        .collect(Collectors.toCollection(TreeSet::new));
  }
}
