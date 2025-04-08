package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.junit.jupiter.api.TestInfo;
import org.opentest4j.AssertionFailedError;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DiplotypeTest;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;

import static org.junit.jupiter.api.Assertions.*;


/**
 * This class wraps {@link Pipeline} to simplify usage in tests.
 *
 * @author Mark Woon
 */
public class PipelineWrapper {
  // controls to support running PipelineTest from SyntheticBatchTest
  private static boolean m_compact = true;
  private static List<PrescribingGuidanceSource> m_sources = PrescribingGuidanceSource.listValues();

  private final Env m_env;
  private final Path m_outputPath;
  private final TestVcfBuilder m_vcfBuilder;
  private final boolean m_findCombinations;
  private final boolean m_callCyp2d6;
  private final boolean m_topCandidatesOnly;
  private boolean m_compactReport = m_compact;
  private boolean m_deleteIntermediateFiles = m_compact;
  private ReportContext m_reportContext;


  static void setCompact(boolean compact) {
    m_compact = compact;
  }

  static void setSources(List<PrescribingGuidanceSource> sources) {
    m_sources = sources;
  }


  public PipelineWrapper(TestInfo testInfo, boolean allMatches) throws IOException, ReportableException {
    this(testInfo, false, allMatches, false);
  }

  public PipelineWrapper(TestInfo testInfo, boolean findCombinations, boolean allMatches, boolean callCyp2d6)
      throws IOException, ReportableException {
    this(testInfo, null, findCombinations, allMatches, callCyp2d6);
  }

  public PipelineWrapper(TestInfo testInfo, @Nullable String name, boolean findCombinations, boolean allMatches,
      boolean callCyp2d6) throws IOException, ReportableException {
    Preconditions.checkNotNull(testInfo);

    m_env = new Env();
    m_outputPath = TestUtils.getTestOutputDir(testInfo, false);
    if (!Files.isDirectory(m_outputPath)) {
      Files.createDirectories(m_outputPath);
    }
    if (name == null) {
      m_vcfBuilder = new TestVcfBuilder(testInfo).saveFile();
    } else {
      m_vcfBuilder = new TestVcfBuilder(testInfo, name).saveFile();
    }
    m_findCombinations = findCombinations;
    m_callCyp2d6 = callCyp2d6;
    m_topCandidatesOnly = !allMatches;
  }

  PipelineWrapper extendedReport() {
    m_compactReport = false;
    return this;
  }

  public PipelineWrapper saveIntermediateFiles() {
    m_deleteIntermediateFiles = false;
    return this;
  }

  ReportContext getContext() {
    return m_reportContext;
  }

  public TestVcfBuilder getVcfBuilder() {
    return m_vcfBuilder;
  }


  public @Nullable Path executeWithOutsideCalls(Path... outsideCallPath) throws Exception {
    if (outsideCallPath == null || outsideCallPath.length == 0) {
      return execute();
    }
    return execute(null, ImmutableList.copyOf(outsideCallPath), null, false);
  }

  /**
   * Execute the pipeline with the specified VCF file (this file will be copied to the output location).
   *
   * @return path to actual VCF used
   */
  public Path executeWithVcf(Path vcfFile) throws Exception {
    return execute(vcfFile, null, null, false);
  }

  public @Nullable Path execute() throws Exception {
    return execute(null, null, null, false);
  }


  public @Nullable Path execute(@Nullable Path vcfFile, @Nullable List<Path> outsideCallPaths,
      Path sampleMetadataFile, boolean allowNoData) throws Exception {
    VcfFile vcfFileObj = null;
    boolean runMatcher = false;
    if (vcfFile != null) {
      runMatcher = true;
      Path copy = m_outputPath.resolve(vcfFile.getFileName());
      Files.copy(vcfFile, copy, StandardCopyOption.REPLACE_EXISTING);
      vcfFile = copy;
      vcfFileObj = new VcfFile(vcfFile, false);
    } else if (m_vcfBuilder.hasData() || allowNoData) {
      runMatcher = true;
      vcfFile = m_vcfBuilder.generate();
      vcfFileObj = new VcfFile(vcfFile, false);
    }
    Pipeline pcat = new Pipeline(m_env,
        runMatcher, vcfFileObj, null, true,
        m_topCandidatesOnly, m_callCyp2d6, m_findCombinations, true,
        true, null, outsideCallPaths,
        true, null, null, m_sources, m_compactReport, true, true, true,
        m_outputPath, null, m_deleteIntermediateFiles,
        Pipeline.Mode.TEST, null, false, sampleMetadataFile
    );
    pcat.call();
    m_reportContext = pcat.getReportContext();
    return vcfFile;
  }


  public Env getEnv() {
    return m_env;
  }

  public Path getOutputDir() {
    return m_outputPath;
  }


  private List<String> stripHomozygousNotes(List<String> calls) {
    return calls.stream()
        .map(d -> {
          if (d.endsWith(" (homozygous)")) {
            return d.substring(0, d.length() - 13);
          }
          return d;
        })
        .toList();
  }

  /**
   * Tests the original diplotype from the source (i.e. matcher or outside call).
   */
  void testSourceDiplotypes(DataSource source, String gene, List<String> diplotypes) {
    testSourceDiplotypes(source, gene, diplotypes, null);
  }

  /**
   * Tests the original diplotype from the source (i.e. matcher or outside call).
   *
   * @param visibleAlleleDiplotypes - If non-null, this is compared against diplotype labels and
   * {@code diplotypes} will be compared against allele names. Normally used to test CPIC-style names or outside calls
   * with missing diplotype.
   */
  void testSourceDiplotypes(DataSource source, String gene, List<String> diplotypes,
      @Nullable List<String> visibleAlleleDiplotypes) {
    GeneReport geneReport = getContext().getGeneReport(source, gene);
    assertNotNull(geneReport);

    if (visibleAlleleDiplotypes != null) {
      testDiplotypeLabels(geneReport, visibleAlleleDiplotypes);
      testDiplotypeNames(geneReport, diplotypes);
    } else {
      testDiplotypeLabels(geneReport, diplotypes);
    }

  }

  private void testDiplotypeLabels(GeneReport geneReport, List<String> diplotypes) {
    List<String> actualDiplotypes = geneReport.getSourceDiplotypes().stream()
        .map(Diplotype::getLabel)
        .toList();
    List<String> expectedDiplotypes = stripHomozygousNotes(diplotypes);
    try {
      assertEquals(expectedDiplotypes, actualDiplotypes);
    } catch (AssertionError ex) {
      System.out.println(printDiagnostic(geneReport));
      throw ex;
    }
  }

  private void testDiplotypeNames(GeneReport geneReport, List<String> diplotypes) {
    List<String> actualDiplotypes = geneReport.getSourceDiplotypes().stream()
        .map(d -> {
          if (d.getAllele1() == null) {
            return Haplotype.UNKNOWN + TextConstants.GENOTYPE_DELIMITER + Haplotype.UNKNOWN;
          }
          List<String> haps = new ArrayList<>();
          haps.add(d.getAllele1().getName());
          if (d.getAllele2() != null) {
            haps.add(d.getAllele2().getName());
          }
          haps.sort(HaplotypeNameComparator.getComparator());
          return String.join(TextConstants.GENOTYPE_DELIMITER, haps);
        })
        .toList();
    List<String> expectedDiplotypes = stripHomozygousNotes(diplotypes);
    try {
      assertEquals(expectedDiplotypes, actualDiplotypes);
    } catch (AssertionError ex) {
      System.out.println(printDiagnostic(geneReport));
      throw ex;
    }
  }


  /**
   * Test the "print" calls for a gene that will display in the final report or in the phenotyper. This will check
   * that the call count matches and then check that each call is present (can be 1 or more).
   *
   * @param gene the gene to get diplotypes for
   * @param calls the expected display of the calls, 1 or more
   */
  void testPrintCpicCalls(String gene, String... calls) {
    testPrintCalls(DataSource.CPIC, gene, calls);
  }

  void testPrintCalls(DataSource source, String gene, String... calls) {
    testPrintCalls(source, gene, Arrays.asList(calls));
  }

  void testPrintCalls(DataSource source, String gene, List<String> calls) {
    GeneReport geneReport = getContext().getGeneReport(source, gene);
    assertNotNull(geneReport, "Missing report for " + gene);
    List<String> actualCalls = ReportHelpers.amdGeneCalls(geneReport);
    List<String> expectedCalls = stripHomozygousNotes(calls);
    try {
      assertEquals(expectedCalls, actualCalls, "Mismatched calls for " + gene);
    } catch (AssertionError ex) {
      System.out.println(printDiagnostic(geneReport));
      throw ex;
    }
  }


  /**
   * Test to make sure the given phenotype is in the collection of phenotypes that come from the source diplotype
   * collection.
   *
   * @param source The source for the diplotype
   * @param gene The gene symbol
   * @param phenotype The phenotype string to check for
   */
  void testSourcePhenotype(DataSource source, String gene, String phenotype) {
    GeneReport geneReport = getContext().getGeneReport(source, gene);
    assertNotNull(geneReport);
    Set<String> sourcePhenotypes = geneReport.getSourceDiplotypes().stream()
        .flatMap(d -> d.getPhenotypes().stream())
        .collect(Collectors.toSet());
    assertTrue(sourcePhenotypes.contains(phenotype), sourcePhenotypes + " does not contain " + phenotype);
  }

  void testRecommendedPhenotype(DataSource source, String gene, String phenotype) {
    GeneReport geneReport = getContext().getGeneReport(source, gene);
    assertNotNull(geneReport);
    SortedSet<String> recPhenotypes = geneReport.getRecommendationDiplotypes().stream()
        .flatMap(d -> d.getPhenotypes().stream())
        .collect(Collectors.toCollection(TreeSet::new));
    assertEquals(1, recPhenotypes.size());
    assertEquals(phenotype, recPhenotypes.first());
  }

  void testRecommendedDiplotypes(String gene, String... haplotypes) {
    testRecommendedDiplotypes(DataSource.CPIC, gene, Arrays.asList(haplotypes));
  }

  /**
   * Test the diplotype that will be used for looking up the recommendation. This will mostly match what's printed in
   * displays but will differ for particular genes.
   *
   * @param gene the gene to get diplotypes for
   * @param haplotypes the expected haplotypes names used for calling, specifying one will assume homozygous,
   * otherwise specify two haplotype names
   */
  void testRecommendedDiplotypes(DataSource source, String gene, List<String> haplotypes) {
    Preconditions.checkArgument(!haplotypes.isEmpty() && haplotypes.size() <= 2,
        "Can only test on 1 or 2 haplotypes, got " + haplotypes.size() + ": " +
            String.join(", ", haplotypes));

    GeneReport geneReport = getContext().getGeneReport(source, gene);
    assertNotNull(geneReport);
    if (haplotypes.size() == 1 && haplotypes.get(0).equals("Unknown/Unknown")) {
      assertFalse(geneReport.isReportable());
      return;
    }
    assertTrue(geneReport.isReportable(), "Not reportable: " + geneReport.getRecommendationDiplotypes());

    Map<String, Integer> lookup = new HashMap<>();
    haplotypes = stripHomozygousNotes(haplotypes);
    if (haplotypes.size() == 2) {
      if (haplotypes.get(0).equals(haplotypes.get(1))) {
        lookup.put(haplotypes.get(0), 2);
      } else {
        lookup.put(haplotypes.get(0), 1);
        lookup.put(haplotypes.get(1), 1);
      }
    } else {
      lookup.put(haplotypes.get(0), 1);
    }

    assertTrue(geneReport.getRecommendationDiplotypes().stream()
            .anyMatch(d -> DiplotypeTest.computeLookupMap(d).equals(lookup)),
        "EXPECTED: " + lookup + " -- ACTUAL: " +
            geneReport.getRecommendationDiplotypes().stream().map(DiplotypeTest::computeLookupMap).toList());
  }

  void testDpydLookup(List<String> haplotypes) {

    GeneReport geneReport = getContext().getGeneReport(DataSource.CPIC, "DPYD");
    assertNotNull(geneReport);
    assertTrue(geneReport.isReportable(), "Not reportable: " + geneReport.getRecommendationDiplotypes());

    System.out.println(geneReport.getRecommendationDiplotypes());
    geneReport.getRecommendationDiplotypes()
        .forEach(d -> System.out.println(d + " -- " + DiplotypeTest.computeLookupMap(d)));

  }

  void testLookupByActivity(DataSource source, String gene, String activityScore) {
    GeneReport geneReport = getContext().getGeneReport(source, gene);
    assertNotNull(geneReport);
    assertTrue(geneReport.isReportable());
    String foundScores = geneReport.getRecommendationDiplotypes().stream()
        .map(Diplotype::getActivityScore)
        .collect(Collectors.joining("; "));
    assertTrue(geneReport.getRecommendationDiplotypes().stream()
            .allMatch(d -> printLookupKeys(d).equals(activityScore)),
        "Activity score mismatch, expected " + activityScore + " but got " + foundScores);
  }

  /**
   * Check to see if all the given genes have been called by the matcher
   */
  void testCalledByMatcher(String... genes) {
    assertTrue(genes != null && genes.length > 0);
    Arrays.stream(genes)
        .forEach(g -> {
          assertTrue(getContext().getGeneReports(g).stream().allMatch(GeneReport::isCalled) &&
                  getContext().getGeneReports(g).stream().noneMatch(GeneReport::isOutsideCall),
              g + " is not called");
        });
  }

  void testReportable(String... genes) {
    assertTrue(genes != null && genes.length > 0);
    Arrays.stream(genes)
        .forEach(g -> {
          assertTrue(getContext().getGeneReports(g).stream().anyMatch(GeneReport::isReportable),
              g + " is not reportable");
        });
  }

  /**
   * Check to see if none of the given genes have been called by the matcher
   */
  void testNotCalledByMatcher(String... genes) {
    Preconditions.checkArgument(genes != null && genes.length > 0);
    for (String g : genes) {
      try {
        assertTrue(getContext().getGeneReports(g).stream().allMatch(GeneReport::isOutsideCall) ||
                getContext().getGeneReports(g).stream().noneMatch(GeneReport::isCalled),
            g + " is called");
      } catch (AssertionFailedError ex) {
        for (GeneReport gr : getContext().getGeneReports(g)) {
          System.out.println(printDiagnostic(gr));
        }
        throw ex;
      }
    }
  }

  /**
   * Check to see if there is a matching recommendation for the given drug name.
   *
   * @param drugName a drug name that has recommendations
   * @param expectedCount the number of matching recommendations you expect
   */
  void testMatchedAnnotations(String drugName, int expectedCount) {
    List<DrugReport> drugReports = getContext().getDrugReports(drugName);
    int numMatched = drugReports.stream()
        .mapToInt(DrugReport::getMatchedAnnotationCount)
        .sum();
    assertEquals(expectedCount, numMatched,
        drugName + " has " + numMatched + " matching recommendation(s) instead of " + expectedCount);
  }

  void testMatchedAnnotations(String drugName, PrescribingGuidanceSource source, int expectedCount) {
    DrugReport drugReport = getContext().getDrugReport(source, drugName);
    assertNotNull(drugReport);
    assertEquals(expectedCount, drugReport.getMatchedAnnotationCount(),
        drugName + " has " + drugReport.getMatchedAnnotationCount() + " matching " + source +
            " recommendation(s) instead of " + expectedCount);
  }

  void testAnyMatchFromSource(String drugName, PrescribingGuidanceSource source) {
    DrugReport drugReport = getContext().getDrugReport(source, drugName);
    assertNotNull(drugReport);
    assertTrue(drugReport.getGuidelines().stream().anyMatch((g) -> g.getSource() == source && g.isMatched()),
        drugName + " does not have matching recommendation from " + source);
  }

  void testNoMatchFromSource(String drugName, PrescribingGuidanceSource source) {
    DrugReport drugReport = getContext().getDrugReport(source, drugName);
    if (drugReport != null) {
      assertTrue(drugReport.getGuidelines().stream().noneMatch(r -> r.getSource() == source && r.isMatched()),
          drugName + " has a matching recommendation from " + source + " and expected none");
    }
  }

  void testMessageCountForDrug(PrescribingGuidanceSource source, String drugName, int messageCount) {
    DrugReport drugReport = getContext().getDrugReport(source, drugName);
    assertNotNull(drugReport);
    assertEquals(messageCount, drugReport.getMessages().size(),
        drugName + " expected " + messageCount + " messages and got " + drugReport.getMessages());
  }

  void testMessageCountForGene(DataSource source, String geneName, int messageCount) {
    GeneReport report = getContext().getGeneReport(source, geneName);
    assertNotNull(report);
    assertEquals(messageCount, report.getMessages().size(),
        geneName + " expected " + messageCount + " messages and got " + report.getMessages());
  }

  void testGeneHasMessage(DataSource source, String geneName, String... msgNames) {
    GeneReport report = getContext().getGeneReport(source, geneName);
    assertNotNull(report);
    assertTrue(report.getMessages().size() >= msgNames.length);
    List<String> names = report.getMessages().stream()
        .map(MessageAnnotation::getName)
        .toList();
    for (String msgName : msgNames) {
      assertTrue(names.contains(msgName),
          geneName + " is missing expected message with name \"" + msgName + "\" (got " + names + ")");
    }
  }


  private String printDiagnostic(GeneReport geneReport) {
    return String.format("""
        
        Matcher:              %s
        Reporter:             %s
        Print (displayCalls): %s
        """,
        geneReport.getSourceDiplotypes().toString(),
        geneReport.getRecommendationDiplotypes().toString(),
        ReportHelpers.amdGeneCalls(geneReport)
    );
  }

  public String printLookupKeys(Diplotype diplotype) {
    return String.join(";", diplotype.getLookupKeys());
  }
}
