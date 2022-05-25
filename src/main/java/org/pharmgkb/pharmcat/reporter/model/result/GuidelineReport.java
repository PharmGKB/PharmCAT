package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;


public class GuidelineReport {

  @Expose
  @SerializedName("name")
  private String name;
  @Expose
  @SerializedName("id")
  private String id;
  @Expose
  @SerializedName("source")
  private DataSource source;
  @Expose
  @SerializedName("version")
  private String version;
  @Expose
  @SerializedName("url")
  private String url;
  @Expose
  @SerializedName("annotationGroups")
  private List<AnnotationGroup> annotationGroups = new ArrayList<>();
  private transient final SortedSet<GeneReport> relatedGeneReports = new TreeSet<>();

  public GuidelineReport(Drug cpicDrug) {
    setVersion(cpicDrug.getCpicVersion());
    setId(cpicDrug.getDrugId());
    setName(cpicDrug.getGuidelineName());
    setSource(DataSource.CPIC);
    setUrl(cpicDrug.getUrl());
  }

  public GuidelineReport(GuidelinePackage guidelinePackage) {
    setVersion(String.valueOf(guidelinePackage.getVersion()));
    setId(guidelinePackage.getGuideline().getId());
    setName(guidelinePackage.getGuideline().getName());
    setSource(DataSource.DPWG);
    setUrl(guidelinePackage.getGuideline().getUrl());
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }

  public String getId() {
    return id;
  }

  public void setId(String id) {
    this.id = id;
  }

  public DataSource getSource() {
    return source;
  }

  public void setSource(DataSource source) {
    this.source = source;
  }

  public String getVersion() {
    return version;
  }

  public void setVersion(String version) {
    this.version = version;
  }

  public String getUrl() {
    return url;
  }

  public void setUrl(String url) {
    this.url = url;
  }

  public List<AnnotationGroup> getAnnotationGroups() {
    return annotationGroups;
  }

  public void setAnnotationGroups(List<AnnotationGroup> annotationGroups) {
    this.annotationGroups = annotationGroups;
  }

  public void addAnnotationGroup(AnnotationGroup annotationGroup) {
    annotationGroups.add(annotationGroup);
  }

  public boolean isMatched() {
    return annotationGroups.size() > 0;
  }

  public Set<String> getRelatedGenes() {
    return getRelatedGeneReports().stream().map(GeneReport::getGene).collect(Collectors.toSet());
  }

  public SortedSet<GeneReport> getRelatedGeneReports() {
    return relatedGeneReports;
  }

  public void addRelatedGeneReport(GeneReport geneReport) {
    if (geneReport != null) {
      relatedGeneReports.add(geneReport);
    }
  }

  public String getUncalledGenes() {
    return relatedGeneReports.stream()
        .filter(g -> !g.isCalled())
        .map(GeneReport::getGene)
        .collect(Collectors.joining(", "));
  }

  public boolean isUncallable() {
    return relatedGeneReports.stream()
        .noneMatch(GeneReport::isCalled);
  }
}
