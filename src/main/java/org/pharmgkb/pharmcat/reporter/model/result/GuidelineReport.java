package org.pharmgkb.pharmcat.reporter.model.result;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import com.google.common.collect.HashMultimap;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.cpic.Recommendation;
import org.pharmgkb.pharmcat.reporter.model.pgkb.Group;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;


public class GuidelineReport {
  @Expose
  @SerializedName("id")
  private String m_id;
  @Expose
  @SerializedName("name")
  private String m_name;

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

  private transient final SortedSet<String> m_genes = new TreeSet<>();
  private transient final SortedSet<GeneReport> m_geneReports = new TreeSet<>();
  private transient List<Genotype> m_possibleGenotypes;


  /**
   * Private constructor for GSON.
   */
  private GuidelineReport() {
  }

  public GuidelineReport(Drug cpicDrug, ReportContext reportContext) {
    m_id = cpicDrug.getDrugId();
    m_name = cpicDrug.getGuidelineName();
    m_source = DataSource.CPIC;
    m_version = cpicDrug.getCpicVersion();
    m_url = cpicDrug.getUrl();
    initializeGenes(cpicDrug.getGenes(), reportContext);
    matchAnnotations(cpicDrug);
  }

  public GuidelineReport(GuidelinePackage guidelinePackage, ReportContext reportContext) {
    m_id = guidelinePackage.getGuideline().getId();
    m_name = guidelinePackage.getGuideline().getName();
    m_source = DataSource.DPWG;
    m_version = guidelinePackage.getVersion();
    m_url = guidelinePackage.getGuideline().getUrl();
    initializeGenes(guidelinePackage.getGenes(), reportContext);
    matchAnnotations(guidelinePackage);
  }

  private void initializeGenes(Collection<String> genes, ReportContext reportContext) {
    // link guideline report to gene report
    for (String geneSymbol : genes) {
      GeneReport geneReport = reportContext.getGeneReport(m_source, geneSymbol);
      if (geneReport == null) {
        continue;
      }
      m_genes.add(geneSymbol);
      m_geneReports.add(geneReport);
    }
    m_possibleGenotypes = Genotype.makeGenotypes(m_geneReports);
  }



  public String getId() {
    return m_id;
  }

  public String getName() {
    return m_name;
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


  public SortedSet<String> getGenes() {
    return m_genes;
  }

  public SortedSet<GeneReport> getGeneReports() {
    return m_geneReports;
  }

  public List<String> getUncalledGenes() {
    return m_geneReports.stream()
        .filter(g -> !g.isCalled())
        .map(GeneReport::getGeneDisplay)
        .toList();
  }

  public boolean isUncallable() {
    return m_geneReports.stream()
        .noneMatch(GeneReport::isCalled);
  }


  private void matchAnnotations(Drug cpicDrug) {
    if (cpicDrug.getDrugName().equals("warfarin")) {
      AnnotationReport annGroup = AnnotationReport.forWarfarin(m_possibleGenotypes);
      m_annotationReports.add(annGroup);
    } else if (cpicDrug.getRecommendations() != null) {
      int gx = 0;
      Map<Recommendation, AnnotationReport> recommendationMap = new HashMap<>();
      for (Genotype genotype : m_possibleGenotypes) {
        gx += 1;
        for (Recommendation rec : cpicDrug.getRecommendations()) {
          if (rec.matchesGenotype(genotype)) {
            AnnotationReport annotationReport = recommendationMap.get(rec);
            if (annotationReport == null) {
              String id = "cpic-" + cpicDrug.getDrugName() + "-" + (recommendationMap.size() + 1) + "-" + gx;
              annotationReport = new AnnotationReport(rec, id);
              recommendationMap.put(rec, annotationReport);
            }
            annotationReport.addGenotype(genotype);
          }
        }
      }
      recommendationMap.values().forEach((ar) -> {
        ar.checkDiplotypes();
        m_annotationReports.add(ar);
      });
    }
   }


  private void matchAnnotations(GuidelinePackage guidelinePackage) {
    HashMultimap<Group, Genotype> matchedGenotypes = HashMultimap.create();
    for (Genotype genotype : m_possibleGenotypes) {
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
      AnnotationReport annotationReport = new AnnotationReport(group, geneSymbol, id);
      matchedGenotypes.get(group).forEach((genotype) -> {
        // TODO(markwoon): check with Ryan why we're doing this
        guidelinePackage.applyFunctions(genotype);
        annotationReport.addGenotype(genotype);
      });
      annotationReport.checkDiplotypes();
      m_annotationReports.add(annotationReport);
    }
  }
}
