package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;


public class GuidelineReport {

  private String name;
  private String id;
  private DataSource source;
  private String version;
  private String url;
  private List<AnnotationGroup> annotationGroups = new ArrayList<>();
  private final SortedSet<String> relatedGenes = new TreeSet<>();
  private final SortedSet<String> uncalledGenes = new TreeSet<>();

  public GuidelineReport(Drug cpicDrug) {
    setVersion(cpicDrug.getCpicVersion());
    setId(cpicDrug.getDrugId());
    setName(cpicDrug.getGuidelineName());
    setSource(DataSource.CPIC);
    setUrl(cpicDrug.getUrl());
    relatedGenes.addAll(cpicDrug.getGenes());
  }

  public GuidelineReport(GuidelinePackage guidelinePackage) {
    setVersion(String.valueOf(guidelinePackage.getVersion()));
    setId(guidelinePackage.getGuideline().getId());
    setName(guidelinePackage.getGuideline().getName());
    setSource(DataSource.DPWG);
    setUrl(guidelinePackage.getGuideline().getUrl());
    relatedGenes.addAll(guidelinePackage.getGenes());
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

  public SortedSet<String> getRelatedGenes() {
    return relatedGenes;
  }

  public boolean isReportable() {
    return this.uncalledGenes.isEmpty();
  }

  public SortedSet<String> getUncalledGenes() {
    return this.uncalledGenes;
  }

  public void addUncalledGene(String uncalledGene) {
    this.uncalledGenes.add(uncalledGene);
  }
}
