package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.collect.HashMultimap;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.pgkb.Group;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;


public class GuidelineReport {
  @Expose
  @SerializedName("name")
  private String m_name;
  @Expose
  @SerializedName("id")
  private String m_id;
  @Expose
  @SerializedName("source")
  private DataSource m_source;
  @Expose
  @SerializedName("version")
  private String m_version;
  @Expose
  @SerializedName("url")
  private String m_url;
  @Expose
  @SerializedName("annotations")
  private List<AnnotationReport> m_annotationReports = new ArrayList<>();

  private transient final SortedSet<GeneReport> m_relatedGeneReports = new TreeSet<>();


  public GuidelineReport(Drug cpicDrug) {
    m_version = cpicDrug.getCpicVersion();
    m_id = cpicDrug.getDrugId();
    m_name = cpicDrug.getGuidelineName();
    m_source = DataSource.CPIC;
    m_url = cpicDrug.getUrl();
  }

  public GuidelineReport(GuidelinePackage guidelinePackage) {
    m_version = String.valueOf(guidelinePackage.getVersion());
    m_id = guidelinePackage.getGuideline().getId();
    m_name = guidelinePackage.getGuideline().getName();
    m_source = DataSource.DPWG;
    m_url = guidelinePackage.getGuideline().getUrl();
  }


  public String getName() {
    return m_name;
  }

  public String getId() {
    return m_id;
  }

  public DataSource getSource() {
    return m_source;
  }

  public String getVersion() {
    return m_version;
  }

  public String getUrl() {
    return m_url;
  }


  public List<AnnotationReport> getAnnotations() {
    return m_annotationReports;
  }

  public boolean isMatched() {
    return m_annotationReports.size() > 0;
  }


  public Set<String> getRelatedGenes() {
    return m_relatedGeneReports.stream().map(GeneReport::getGene).collect(Collectors.toSet());
  }

  public SortedSet<GeneReport> getRelatedGeneReports() {
    return m_relatedGeneReports;
  }

  public void addRelatedGeneReport(GeneReport geneReport) {
    if (geneReport != null) {
      m_relatedGeneReports.add(geneReport);
    }
  }

  public List<String> getUncalledGenes() {
    return m_relatedGeneReports.stream()
        .filter(g -> !g.isCalled())
        .map(GeneReport::getGeneDisplay)
        .toList();
  }

  public boolean isUncallable() {
    return m_relatedGeneReports.stream()
        .noneMatch(GeneReport::isCalled);
  }



  public void matchAnnotationsToGenotype(List<Genotype> genotypes, Drug cpicDrug) {
    if (cpicDrug.getDrugName().equals("warfarin")) {
      AnnotationReport annGroup = AnnotationReport.forWarfarin(genotypes);
      m_annotationReports.add(annGroup);
    } else if (cpicDrug.getRecommendations() != null) {
      int gx = 0;
      for (Genotype genotype : genotypes) {
        gx += 1;
        int rx = 0;
        for (Recommendation rec : cpicDrug.getRecommendations()) {
          if (rec.matchesGenotype(genotype)) {
            rx += 1;
            String id = "cpic-" + cpicDrug.getDrugName() + "-" + rx + "-" + gx;
            AnnotationReport annGroup = new AnnotationReport(rec, id);
            annGroup.addGenotype(genotype);
            m_annotationReports.add(annGroup);
          }
        }
      }
    }
   }


  public void matchAnnotationsToGenotype(List<Genotype> genotypes, GuidelinePackage guidelinePackage) {
    HashMultimap<Group, Genotype> matchedGenotypes = HashMultimap.create();
    for (Genotype genotype : genotypes) {
      for (Diplotype diplotype : genotype.getDiplotypes()) {
        if (diplotype.isPhenotypeOnly() || diplotype.isAllelePresenceType()) {
          guidelinePackage.getGroups().stream()
              .filter(group -> diplotype.getPhenotypes().stream().anyMatch(p -> group.getName().equalsIgnoreCase(p)))
              .forEach(group -> matchedGenotypes.put(group, genotype));
        } else if (!diplotype.isUnknownAlleles()) {
          Set<String> functionKeys = guidelinePackage.getGuideline().getFunctionKeysForDiplotype(diplotype);
          for (String functionKey : functionKeys) {
            guidelinePackage.getGroups().stream()
                .filter(group -> group.matchesKey(functionKey))
                .forEach(group -> matchedGenotypes.put(group, genotype));
          }
        }
      }
    }
    int rx = 0;
    for (Group group : matchedGenotypes.keys()) {
      rx += 1;
      String geneSymbol = guidelinePackage.getGenes().iterator().next();
      String id = "dpwg-" + guidelinePackage.getGuideline().getId() + "-" + rx;
      AnnotationReport annGroup = new AnnotationReport(group, geneSymbol, id);
      matchedGenotypes.get(group).forEach((genotype) -> {
        guidelinePackage.applyFunctions(genotype);
        annGroup.addGenotype(genotype);
      });
      m_annotationReports.add(annGroup);
    }
  }
}
