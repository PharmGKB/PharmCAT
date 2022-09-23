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
import org.pharmgkb.pharmcat.reporter.model.pgkb.Group;
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
  private List<AnnotationReport> m_annotationReports = new ArrayList<>();

  private transient final SortedSet<GeneReport> relatedGeneReports = new TreeSet<>();


  public GuidelineReport(Drug cpicDrug) {
    version = cpicDrug.getCpicVersion();
    id = cpicDrug.getDrugId();
    name = cpicDrug.getGuidelineName();
    source = DataSource.CPIC;
    url = cpicDrug.getUrl();
  }

  public GuidelineReport(GuidelinePackage guidelinePackage) {
    version = String.valueOf(guidelinePackage.getVersion());
    id = guidelinePackage.getGuideline().getId();
    name = guidelinePackage.getGuideline().getName();
    source = DataSource.DPWG;
    url = guidelinePackage.getGuideline().getUrl();
  }


  public String getName() {
    return name;
  }

  public String getId() {
    return id;
  }

  public DataSource getSource() {
    return source;
  }

  public String getVersion() {
    return version;
  }

  public String getUrl() {
    return url;
  }


  public List<AnnotationReport> getAnnotations() {
    return m_annotationReports;
  }

  public void addAnnotation(AnnotationReport annotationReport) {
    m_annotationReports.add(annotationReport);
  }

  public boolean isMatched() {
    return m_annotationReports.size() > 0;
  }


  public Set<String> getRelatedGenes() {
    return relatedGeneReports.stream().map(GeneReport::getGene).collect(Collectors.toSet());
  }

  public SortedSet<GeneReport> getRelatedGeneReports() {
    return relatedGeneReports;
  }

  public void addRelatedGeneReport(GeneReport geneReport) {
    if (geneReport != null) {
      relatedGeneReports.add(geneReport);
    }
  }

  public List<String> getUncalledGenes() {
    return relatedGeneReports.stream()
        .filter(g -> !g.isCalled())
        .map(GeneReport::getGeneDisplay)
        .toList();
  }

  public boolean isUncallable() {
    return relatedGeneReports.stream()
        .noneMatch(GeneReport::isCalled);
  }



  public void matchAnnotationsToGenotype(List<Genotype> genotypes, Drug cpicDrug) {
    if (cpicDrug.getDrugName().equals("warfarin")) {
      AnnotationReport annGroup = AnnotationReport.forWarfarin(genotypes);
      addAnnotation(annGroup);
    } else if (cpicDrug.getRecommendations() != null) {
      for (Genotype genotype : genotypes) {
        cpicDrug.getRecommendations().stream()
            .filter((r) -> r.matchesGenotype(genotype))
            .forEach((r) -> {
              AnnotationReport annGroup = new AnnotationReport(r);
              annGroup.addGenotype(genotype);
              addAnnotation(annGroup);
            });
      }
    }
  }


  public void matchAnnotationsToGenotype(List<Genotype> genotypes, GuidelinePackage guidelinePackage) {
    Set<Group> matchedGroups = new TreeSet<>();
    for (Genotype genotype : genotypes) {
      for (Diplotype diplotype : genotype.getDiplotypes()) {
        if (diplotype.isPhenotypeOnly() || diplotype.isAllelePresenceType()) {
          guidelinePackage.getGroups().stream()
              .filter(group -> diplotype.getPhenotypes().stream().anyMatch(p -> group.getName().equalsIgnoreCase(p)))
              .forEach(g -> {
                g.addMatchingDiplotype(diplotype);
                g.addMatchingGenotype(genotype);
                matchedGroups.add(g);
              });
        } else if (!diplotype.isUnknownAlleles()) {
          Set<String> functionKeys = guidelinePackage.getGuideline().getFunctionKeysForDiplotype(diplotype);
          for (String functionKey : functionKeys) {
            guidelinePackage.getGroups().stream()
                .filter(group -> group.matchesKey(functionKey))
                .forEach(g -> {
                  g.addMatchingFunctionKey(functionKey);
                  g.addMatchingDiplotype(diplotype);
                  g.addMatchingGenotype(genotype);
                  matchedGroups.add(g);
                });
          }
        }
      }
    }
    matchedGroups
        .forEach((group) -> {
          AnnotationReport annGroup = new AnnotationReport(group, guidelinePackage.getGenes().iterator().next());
          group.getMatchingGenotypes().forEach((genotype) -> {
            guidelinePackage.applyFunctions(genotype);
            annGroup.addGenotype(genotype);
          });
          addAnnotation(annGroup);
        });
  }
}
