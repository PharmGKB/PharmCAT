package org.pharmgkb.pharmcat.reporter;

import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.definition.IncidentalFinder;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.AstrolabeCall;
import org.pharmgkb.pharmcat.reporter.model.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.RelatedGene;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReportUgt1a1;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;


/**
 * This class acts as a central context for all data needed to generate the final report.
 *
 * It currently gathers
 * <ul>
 *   <li>{@link GeneCall} objects from the named allele matcher</li>
 *   <li>{@link GuidelineReport} objects from dosing guideline annotations</li>
 *   <li>Allele definitions on a per-gene basis</li>
 * </ul>
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class ReportContext {

  private Map<String,GeneReport> m_geneReports = new TreeMap<>();
  private List<GuidelineReport> m_guidelineReports;
  private PhenotypeMap m_phenotypeMap = new PhenotypeMap();
  private IncidentalFinder m_incidentalFinder = new IncidentalFinder();

  private final Predicate<String> isGeneIncidental = s -> m_geneReports.values().stream()
      .anyMatch(r -> r.getGene().equals(s) && r.isIncidental());

  /**
   * Public constructor. Compiles all the incoming data into useful objects to be held for later reporting
   * @param calls {@link GeneCall} objects from the sample data
   * @param astrolabeCalls {@link AstrolabeCall} objects from the astrolabe data, non-null but can be empty
   * @param guidelinePackages a List of all the guidelines to try to apply
   */
  public ReportContext(List<GeneCall> calls, @Nonnull List<AstrolabeCall> astrolabeCalls, List<GuidelinePackage> guidelinePackages) throws Exception {

    makeGuidelineReports(guidelinePackages);

    makeGeneReports(guidelinePackages);
    compileMatcherData(calls);
    compileAstrolabeData(astrolabeCalls);

    findMatches();

    compileDrugsForGenes();
  }

  /**
   * Applies the given {@link MessageAnnotation} objects to the data in this report.
   * @param messages a List of {@link MessageAnnotation} objects
   */
  public void applyMessage(List<MessageAnnotation> messages) {
    m_geneReports.values().forEach(r -> r.applyMessages(messages));
    m_guidelineReports.forEach(r -> r.applyMessages(messages));
  }

  /**
   * Takes the raw GuidelinePackage objects from PharmGKG and maps them to {@link GuidelineReport} objects that can be
   * used in the reporter
   * @param guidelinePackages a List of PharmGKB {@link GuidelinePackage} objects
   */
  private void makeGuidelineReports(List<GuidelinePackage> guidelinePackages) {
    m_guidelineReports = guidelinePackages.stream().map(GuidelineReport::new).collect(Collectors.toList());
  }

  /**
   * Makes {@link GeneReport} objects for each of the genes found in PharmGKB {@link GuidelinePackage} objects
   * @param guidelinePackages a List of PharmGKB {@link GuidelinePackage} objects
   */
  private void  makeGeneReports(List<GuidelinePackage> guidelinePackages) {
    for (GuidelinePackage guidelinePackage : guidelinePackages) {
      guidelinePackage.getGuideline().getRelatedGenes().stream()
          .map(RelatedGene::getSymbol)
          .distinct()
          .forEach(s -> m_geneReports.put(s, GeneReportFactory.newReport(s)));
    }
  }

  /**
   * Takes {@link GeneCall} data and preps internal data structures for usage. Also prepares exception logic and applies
   * it to the calling data
   */
  private void compileMatcherData(List<GeneCall> calls) throws Exception {
    for (GeneCall call : calls) {
      GeneReport geneReport = GeneReportFactory.newReport(call);
      geneReport.setCallData(call);
      m_geneReports.put(call.getGene(), geneReport);

      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(call.getGene(), call.getVariants(), m_phenotypeMap, m_incidentalFinder);
      if (geneReport instanceof GeneReportUgt1a1) {
        GeneReportUgt1a1 ugt1a1Report = (GeneReportUgt1a1)geneReport;
        diplotypeFactory.makeDiplotypes(ugt1a1Report).forEach(geneReport::addDiplotype);
        ugt1a1Report.getDistinctAlleles().stream()
            .map(diplotypeFactory::makeHaplotype)
            .map(Haplotype::getFunction)
            .distinct()
            .forEach(ugt1a1Report::addDisplayFunction);
      }
      else {
        diplotypeFactory.makeDiplotypes(call).forEach(geneReport::addDiplotype);
      }
    }
  }

  /**
   * Takes astrolabe calls, find the GeneReport for each one and then adds astrolabe information to it
   * @param calls astrolabe calls
   */
  private void compileAstrolabeData(List<AstrolabeCall> calls) throws Exception {
    for (AstrolabeCall astrolabeCall : calls) {
      GeneReport geneReport = m_geneReports.get(astrolabeCall.getGene());
      geneReport.setAstrolabeData(astrolabeCall);

      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(astrolabeCall.getGene(), null, m_phenotypeMap, m_incidentalFinder);
      diplotypeFactory.makeDiplotypes(astrolabeCall).forEach(geneReport::addDiplotype);
    }
  }

  /**
   * Adds compiled drug data to each {@link GeneReport} object that's related through a guideline. This gives a useful
   * collection of all drugs that are related to each gene.
   */
  private void compileDrugsForGenes() {
    m_guidelineReports.forEach(r -> {
      for (String gene : r.getRelatedGeneSymbols()) {
        getGeneReport(gene).addRelatedDrugs(r);
      }
    });
  }

  private boolean isCalled(String gene) {
    return m_geneReports.get(gene).isCalled();
  }

  /**
   * Assigns matched guideline groups for all guidelines in this report based on called diplotype functions for each
   * gene and guideline combination.
   */
  private void findMatches() throws Exception {

    for(GuidelineReport guideline : m_guidelineReports) {
      boolean reportable = guideline.getRelatedGeneSymbols().stream()
          .anyMatch(this::isCalled);

      guideline.setReportable(reportable);

      guideline.getRelatedGeneSymbols().stream()
          .filter(g -> !isCalled(g))
          .forEach(guideline::addUncalledGene);

      if (!reportable) {
        continue;
      }

      makeAllReportGenotypes(guideline);
      guideline.setIncidentalResult(
          guideline.getRelatedGeneSymbols().stream()
              .anyMatch(isGeneIncidental));
    }
  }

  /**
   * Makes a set of called genotype function Strings (in the form "GENEA:No Function/No Function;GENEB:Normal Function/Normal Function") for the given
   * {@link GuidelineReport} and then adds them to the given {@link GuidelineReport}.
   *
   * <em>note:</em> this needs to stay in the {@link ReportContext} since it relies on referencing possibly multiple
   * {@link GeneReport} objects.
   */
  private void makeAllReportGenotypes(GuidelineReport guidelineReport) {
    Set<String> results = new TreeSet<>();
    for (String symbol : guidelineReport.getRelatedGeneSymbols()) {
      results = makeCalledGenotypes(guidelineReport, symbol, results);
    }
    results.forEach(guidelineReport::addReportGenotype);
  }

  /**
   * Makes the genotype Strings for the given gene <code>symbol</code> and either creates a new Set if the passed
   * <code>results</code> is empty or adds to the existing <code>results</code> if there are already entries
   *
   * @param guidelineReport a GuidelineReport with gene calls
   * @param symbol the gene to generate calls for
   * @param results the existing call list to add to
   * @return a new Set of gene calls with the calls for the specified gene
   */
  private Set<String> makeCalledGenotypes(GuidelineReport guidelineReport, String symbol, @Nonnull Set<String> results) {
    if (results.size() == 0) {
      return getGeneReport(symbol).getDiplotypeLookupKeys().stream()
          .map(guidelineReport::translateToPhenotype)
          .collect(Collectors.toSet());
    }
    else {
      Set<String> newResults = new TreeSet<>();
      for (String geno1 : results) {
        getGeneReport(symbol).getDiplotypeLookupKeys().stream().map(guidelineReport::translateToPhenotype).forEach(
            geno2 -> {
              Set<String> genos = new TreeSet<>();
              genos.add(geno1);
              genos.add(geno2);
              newResults.add(genos.stream().collect(Collectors.joining(";")));
            });
      }
      return newResults;
    }
  }

  public List<GuidelineReport> getGuidelineReports() {
    return m_guidelineReports;
  }

  public Collection<GeneReport> getGeneReports() {
    return m_geneReports.values();
  }

  @Nonnull
  public GeneReport getGeneReport(@Nonnull String geneSymbol) {
    return m_geneReports.get(geneSymbol);
  }
}
