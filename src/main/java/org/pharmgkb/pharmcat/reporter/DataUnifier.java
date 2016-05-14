package org.pharmgkb.pharmcat.reporter;

import java.lang.invoke.MethodHandles;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeMultimap;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.DosingGuideline;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 *
 * This is the primary class and method for matching data the the three different sources.
 * I hate just about everything about how this was done, but for the sake of a quick hack to get
 * reports up and running it will have to do.
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class DataUnifier {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  private List<GeneCall> m_calls;
  private Multimap<String, String> m_sampleGeneToDiplotypeMap = TreeMultimap.create();
  private Map<String, GeneReport> m_symbolToGeneReportMap = new TreeMap<>();
  private List<GuidelineReport> m_interactionList;
  private Set<String> m_calledGenes;

  /**
   * public constructor
   * @param calls GeneCall objects from the sample data
   * @param guidelines a List of all the guidelines to try to apply
   */
  public DataUnifier(List<GeneCall> calls, List<DosingGuideline> guidelines) throws Exception {
    m_calls = calls;
    m_interactionList = guidelines.stream().map(GuidelineReport::new).collect(Collectors.toList());
    compileGeneData();
    findMatches();
  }

  /**
   * Takes {@link GeneCall} data and preps internal data structures for usage. Also prepares exception logic and applies
   * it to the calling data
   * @throws Exception
   */
  private void compileGeneData() throws Exception {
    ExceptionMatcher exceptionMatcher = new ExceptionMatcher();

    for (GeneCall call : m_calls) {

      // starts a new GeneReport based on data in the GeneCall
      GeneReport geneReport = new GeneReport(call);
      // adds exceptions to the GeneReport
      exceptionMatcher.addExceptions(geneReport);

      m_sampleGeneToDiplotypeMap.putAll(geneReport.getGene(), geneReport.getDips());

      m_symbolToGeneReportMap.put(call.getGene(), geneReport);
    }
    fixCyp2c19();
    m_calledGenes = m_calls.stream().map(GeneCall::getGene).collect(Collectors.toSet());
  }

  /**
   * Substitutes "*4" for "*4A" or "*4B" in the map of called diplotypes
   *
   * This is a temporary fix while the allele definitions for CYP2C19 don't match what's been annotated by CPIC
   *
   */
  private void fixCyp2c19() {
    if (m_sampleGeneToDiplotypeMap.keySet().contains("CYP2C19")) {
      List<String> fixedDiplotypes = m_sampleGeneToDiplotypeMap.get("CYP2C19").stream()
          .map(d -> d.replaceAll("\\*4[AB]", "*4"))
          .collect(Collectors.toList());
      m_sampleGeneToDiplotypeMap.removeAll("CYP2C19");
      m_sampleGeneToDiplotypeMap.putAll("CYP2C19", fixedDiplotypes);
    }
  }

  /**
   *  Call to do the actual matching, this should all be broken out into
   *  independent methods so errors are clearly and atomically identified
   *  and handled.
   *
   *  This is going to need to be rethought through and reconstructed
   */
  private void findMatches() throws Exception {

    for(GuidelineReport guideline : m_interactionList) {
      guideline.setReportable(m_calledGenes);
      if (!guideline.isReportable()) {
        sf_logger.warn("Can't annotate guideline {}, it's missing {}",
            guideline.getName(),
            guideline.getRelatedGeneSymbols().stream().filter(s -> !m_calledGenes.contains(s)).collect(Collectors.joining(",")));
        continue;
      }

      sf_logger.info("Able to use {}", guideline.getName());

      Set<String> calledGenotypesForGuideline = makeAllCalledGenotypes(guideline.getRelatedGeneSymbols());

      for (Group annotationGroup : guideline.getGroups()) {
        calledGenotypesForGuideline.stream()
            .filter(calledGenotype -> annotationGroup.getGenotypes().contains(calledGenotype))
            .forEach(calledGenotype -> {
              guideline.addMatchingGroup(annotationGroup);
              guideline.putMatchedDiplotype(annotationGroup.getId(), calledGenotype);
            });
      }
    }
  }

  /**
   * Makes a set of called genotype Strings for the given collection of genes. This can be used later for matching to
   * annotation group genotypes
   * @param geneSymbols the gene symbols to include in the genotype strings
   * @return a Set of string genotype calls in the form "GENEA:*1/*2;GENEB:*2/*3"
   */
  private Set<String> makeAllCalledGenotypes(Collection<String> geneSymbols) {
    Set<String> results = new TreeSet<>();
    for (String symbol : geneSymbols) {
      results = makeCalledGenotypes(symbol, results);
    }
    return results;
  }

  private Set<String> makeCalledGenotypes(String symbol, Set<String> results) {
    if (results.size() == 0) {
      return Sets.newHashSet(m_sampleGeneToDiplotypeMap.get(symbol));
    }
    else {
      Set<String> newResults = new TreeSet<>();
      for (String geno1 : results) {
        for (String geno2 : m_sampleGeneToDiplotypeMap.get(symbol)) {
          Set<String> genos = new TreeSet<>();
          genos.add(geno1);
          genos.add(geno2);
          newResults.add(genos.stream().collect(Collectors.joining(";")));
        }
      }
      return newResults;
    }
  }

  public List<GuidelineReport> getGuidelineResults() {
    return m_interactionList;
  }

  public Map<String, GeneReport> getSymbolToGeneReportMap() {
    return m_symbolToGeneReportMap;
  }
}
