package org.pharmgkb.pharmcat.reporter;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.AstrolabeCall;
import org.pharmgkb.pharmcat.reporter.model.Group;
import org.pharmgkb.pharmcat.reporter.model.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.PharmcatException;
import org.pharmgkb.pharmcat.reporter.model.RelatedGene;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;


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

  private Set<GeneReport> m_geneReports = new TreeSet<>();
  private List<GuidelineReport> m_guidelineReports;
  private PhenotypeMap m_phenotypeMap = new PhenotypeMap();

  public final Function<String,Stream<String>> mapGeneToDiplotypes = s -> m_geneReports.stream()
      .filter(c -> c.getGene().equals(s))
      .flatMap(c -> c.getDiplotypes().stream().map(e -> e + (c.isAstrolabeCall() ? " (Astrolabe)" : "")))
      .map(d -> s + ":" + d);

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
   * Applies the given {@link PharmcatException} objects to the data in this report.
   * @param exceptions a List of {@link PharmcatException} objects
   */
  public void applyException(List<PharmcatException> exceptions) {
    m_geneReports.forEach(r -> r.applyExceptions(exceptions));
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
  private void makeGeneReports(List<GuidelinePackage> guidelinePackages) {
    for (GuidelinePackage guidelinePackage : guidelinePackages) {
      guidelinePackage.getGuideline().getRelatedGenes().stream()
          .map(RelatedGene::getSymbol)
          .distinct()
          .forEach(s -> m_geneReports.add(new GeneReport(s)));
    }
  }

  /**
   * Takes {@link GeneCall} data and preps internal data structures for usage. Also prepares exception logic and applies
   * it to the calling data
   */
  private void compileMatcherData(List<GeneCall> calls) throws Exception {
    for (GeneCall call : calls) {
      GeneReport geneReport = m_geneReports.stream()
          .filter(r -> r.getGene().equals(call.getGene()))
          .reduce((r1,r2) -> { throw new RuntimeException("Didn't expect more than one report"); })
          .orElseThrow(IllegalStateException::new);
      geneReport.setCallData(call, m_phenotypeMap);
    }
  }

  /**
   * Takes astrolabe calls, find the GeneReport for each one and then adds astrolabe information to it
   * @param calls astrolabe calls
   */
  private void compileAstrolabeData(List<AstrolabeCall> calls) throws Exception {
    for (AstrolabeCall astrolabeCall : calls) {
      GeneReport geneReport = m_geneReports.stream()
          .filter(r -> r.getGene().equals(astrolabeCall.getGene()))
          .reduce((r1,r2) -> { throw new RuntimeException("Didn't expect more than one report"); })
          .orElseThrow(IllegalStateException::new);
      geneReport.setAstrolabeData(astrolabeCall, m_phenotypeMap);
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
    return m_geneReports.stream().filter(r -> r.getGene().equals(gene)).allMatch(GeneReport::isCalled);
  }

  /**
   *  Call to do the actual matching, this should all be broken out into
   *  independent methods so errors are clearly and atomically identified
   *  and handled.
   *
   *  This is going to need to be rethought through and reconstructed
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

      Set<String> calledGenotypesForGuideline = makeAllCalledGenotypes(guideline);

      for (Group annotationGroup : guideline.getGroups()) {
        calledGenotypesForGuideline.stream()
            .filter(calledGenotype -> annotationGroup.getGenePhenotypes().contains(calledGenotype))
            .forEach(calledGenotype -> {
              guideline.addMatchingGroup(annotationGroup);
              guideline.putMatchedDiplotype(annotationGroup.getId(), calledGenotype);
            });
      }
    }
  }

  /**
   * Makes a set of called genotype Strings for the given {@link GuidelineReport}. This can be used later for matching
   * to annotation group genotypes
   * @return a Set of string genotype calls in the form "GENEA:*1/*2;GENEB:*2/*3"
   */
  private Set<String> makeAllCalledGenotypes(GuidelineReport guidelineReport) {
    Set<String> results = new TreeSet<>();
    for (String symbol : guidelineReport.getRelatedGeneSymbols()) {
      results = makeCalledGenotypes(guidelineReport, symbol, results);
    }
    return results;
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

  public Set<GeneReport> getGeneReports() {
    return m_geneReports;
  }

  @Nonnull
  private GeneReport getGeneReport(@Nonnull String geneSymbol) {
    return m_geneReports.stream().filter(r -> r.getGene().equals(geneSymbol))
        .reduce((r1,r2) -> { throw new RuntimeException("Duplicate gene reports found"); })
        .orElseThrow(RuntimeException::new);
  }
}
