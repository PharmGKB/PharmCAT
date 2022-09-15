package org.pharmgkb.pharmcat.reporter.model.pgkb;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;

import static org.pharmgkb.pharmcat.reporter.model.result.Diplotype.DELIMITER;


/**
 * PharmGKB Guideline Annotation Model.
 */
public class DosingGuideline {
  @Expose
  @SerializedName("id")
  private String id;
  @Expose
  @SerializedName("name")
  private String name;
  @Expose
  @SerializedName("relatedChemicals")
  private List<AccessionObject> relatedChemicals = new ArrayList<>();
  @Expose
  @SerializedName("relatedGenes")
  private List<AccessionObject> relatedGenes = new ArrayList<>();
  @Expose
  @SerializedName("source")
  private String source;
  @Expose
  @SerializedName("recommendation")
  private boolean m_recommendation;
  @Expose
  @SerializedName("guidelineGenes")
  private List<GuidelineGene> m_guidelineGenes;
  @Expose
  @SerializedName("@id")
  private String url;
  @Expose
  @SerializedName("version")
  private Integer m_version;


  public String getId() {
    return id;
  }

  public String getName() {
    return name;
  }

  public List<AccessionObject> getRelatedChemicals() {
    return relatedChemicals;
  }

  public List<AccessionObject> getRelatedGenes() {
    return relatedGenes;
  }

  public String getSource() {
    return source;
  }


  public boolean isRecommendation() {
    return m_recommendation;
  }

  public void setRecommendation(boolean recommendation) {
    m_recommendation = recommendation;
  }

  public List<GuidelineGene> getGuidelineGenes() {
    return m_guidelineGenes;
  }

  public void setGuidelineGenes(List<GuidelineGene> guidelineGenes) {
    m_guidelineGenes = guidelineGenes;
  }

  public Optional<GuidelineGene> findGuidelineGeneFor(String geneSymbol) {
    return m_guidelineGenes.stream().filter(g -> g.getGene().getSymbol().equals(geneSymbol)).findFirst();
  }

  public String getUrl() {
    return url;
  }

  public void setUrl(String url) {
    this.url = url;
  }

  public Set<String> getFunctionKeysForDiplotype(Diplotype diplotype) {
    Set<String> functionKeys = new TreeSet<>();
    findGuidelineGeneFor(diplotype.getGene()).ifPresent(guidelineGene -> {
      List<String> functions = new ArrayList<>();
      functions.add(guidelineGene.findFunctionForAllele(diplotype.getAllele1()).orElse(TextConstants.UNKNOWN_FUNCTION));
      functions.add(guidelineGene.findFunctionForAllele(diplotype.getAllele2()).orElse(TextConstants.UNKNOWN_FUNCTION));
      functions.sort(Comparator.naturalOrder());
      functionKeys.add(diplotype.getGene() + ":" + String.join(DELIMITER, functions));
    });
    return functionKeys;
  }

  public Integer getVersion() {
    return m_version;
  }

  protected void applyFunctions(Genotype genotype) {
    for (Diplotype diplotype : genotype.getDiplotypes()) {
      findGuidelineGeneFor(diplotype.getGene()).ifPresent((guidelineGene) -> {
        if (diplotype.getAllele1() != null) {
          guidelineGene.findFunctionForAllele(diplotype.getAllele1()).ifPresent((fnName) -> {
            diplotype.getAllele1().setFunction(fnName);
          });
        }
        if (diplotype.getAllele2() != null) {
          guidelineGene.findFunctionForAllele(diplotype.getAllele2()).ifPresent((fnName) -> {
            diplotype.getAllele2().setFunction(fnName);
          });
        }
      });
    }
  }
}
