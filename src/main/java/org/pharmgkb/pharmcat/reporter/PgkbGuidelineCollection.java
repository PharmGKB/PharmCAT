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
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.pgkb.AccessionObject;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.pgkb.PrescribingGuidanceDataset;
import org.pharmgkb.pharmcat.reporter.model.pgkb.Publication;
import org.pharmgkb.pharmcat.util.DataSerializer;


public class PgkbGuidelineCollection {
  public static final String PRESCRIBING_GUIDANCE_FILE_NAME = "prescribing_guidance.json";
  public static final Path GUIDANCE_DATA_PATH =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/reporter/" + PRESCRIBING_GUIDANCE_FILE_NAME);

  @Expose
  @SerializedName("guidelines")
  private final List<GuidelinePackage> f_guidelinePackages = new ArrayList<>();
  private final SortedSetMultimap<String,GuidelinePackage> f_guidelineMap = TreeMultimap.create(String::compareToIgnoreCase, Comparator.naturalOrder());
  private SortedSet<String> m_genes;
  private final String m_version;


  public PgkbGuidelineCollection() throws IOException {
    this(GUIDANCE_DATA_PATH);
  }

  public PgkbGuidelineCollection(Path guidanceDataPath) throws IOException {
    try (BufferedReader br = Files.newBufferedReader(guidanceDataPath)) {
      PrescribingGuidanceDataset dataset = DataSerializer.GSON.fromJson(br, PrescribingGuidanceDataset.class);

      m_version = dataset.getVersion();

      f_guidelinePackages.addAll(dataset.getGuidelinePackages());

      for (GuidelinePackage guidelinePackage : dataset.getGuidelinePackages()) {
        guidelinePackage.getCitations().forEach(Publication::normalize);
        for (AccessionObject chemical : guidelinePackage.getGuideline().getRelatedChemicals()) {
          f_guidelineMap.put(chemical.getName(), guidelinePackage);
        }
      }
    }
  }

  public List<GuidelinePackage> getGuidelinePackages() {
    return f_guidelinePackages;
  }

  public String getVersion() {
    return m_version;
  }

  public List<GuidelinePackage> findGuidelinePackages(String chemicalName, PrescribingGuidanceSource source) {
    return f_guidelineMap.get(chemicalName).stream()
        .filter(p -> p.isDataSourceType(source))
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

  public Set<GuidelinePackage> getGuidelinesFromSource(PrescribingGuidanceSource source) {
    return f_guidelineMap.values().stream()
        .filter(g -> g.isDataSourceType(source))
        .collect(Collectors.toSet());
  }

  public Set<String> getChemicalsUsedInSource(PrescribingGuidanceSource source) {
    return f_guidelineMap.values().stream()
        .map(GuidelinePackage::getGuideline)
        .filter(g -> g.getSource().equalsIgnoreCase(source.getPgkbSource().getPharmgkbName()) && g.getObjCls().equalsIgnoreCase(source.getPgkbObjectType()))
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


  public void serializeToJson(Path file) throws IOException {
    DataSerializer.serializeToJson(this, file);
  }
}
