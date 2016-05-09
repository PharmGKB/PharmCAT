package org.pharmgkb.pharmcat.reporter;

import java.util.Set;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.reporter.model.CPICinteraction;
import org.pharmgkb.pharmcat.reporter.model.RelatedGene;


/**
 *
 * Core for matching called genes and haplotypes to the list of available drug recommendations
 *
 * @author greytwist
 *
 */
public class DrugRecommendationMatcher {
  private Set<String> m_neededGenes;
  private Set<String> m_definedGeneSymbolSet;

  public DrugRecommendationMatcher(@Nonnull Set<String> availableGenes, @Nonnull CPICinteraction guideline) throws Exception {
    if (guideline == null || availableGenes == null) {
      throw new Exception("Missing data for constructor");
    }

    m_neededGenes = guideline.getRelatedGenes().stream()
        .map(RelatedGene::getSymbol)
        .filter(s -> !availableGenes.contains(s))
        .collect(Collectors.toSet());

    m_definedGeneSymbolSet = guideline.getRelatedGenes().stream()
        .map(RelatedGene::getSymbol)
        .collect(Collectors.toSet());
  }

  public boolean matches() {
    return m_neededGenes.isEmpty();
  }

  /**
   * These are the genes that are needed for this guideline but which were not present in the available gene collection
   * @return a Set of gene symbols
   */
  public Set<String> getNeededGenes() {
    return m_neededGenes;
  }

  /**
   * The symbols of the genes used in this guideline
   * @return
   */
  public Set<String> getDefinedGeneSymbolSet() {
    return m_definedGeneSymbolSet;
  }
}
