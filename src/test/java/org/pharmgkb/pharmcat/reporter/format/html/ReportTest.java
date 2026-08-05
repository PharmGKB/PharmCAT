package org.pharmgkb.pharmcat.reporter.format.html;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.OutsideCallParser;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.util.DataSerializer;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * This is a JUnit test for {@link Report}.
 *
 * @author Mark Woon
 */
class ReportTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  /**
   * Regression test: an inferred, unmatched RYR1 diplotype must be flagged with the generic "Inferred" marker
   * ({@code unmatchedInferred}), not the DPYD-specific one ({@code unmatchedDpydInferred}). {@link Report}
   * previously used {@code isLowestFunctionGene()} (true for both DPYD and RYR1) instead of checking specifically
   * for "DPYD".
   *
   * <p>Real CPIC/FDA guideline data for RYR1 happens to cover every phenotype combination exhaustively, so this
   * scenario can't be reproduced by running an actual drug guideline through the pipeline. Instead, this builds a
   * minimal synthetic guideline package (with no recommendations at all) to force the "unmatched" branch of
   * {@link Report}'s constructor, while using a real inferred RYR1 diplotype (via an outside call combining more
   * than 2 haplotypes, the same mechanism used for real RYR1 combination calls).
   */
  @Test
  void unmatchedInferredRyr1NotMislabeledAsDpyd() throws Exception {
    Phenotyper phenotyper = new Phenotyper(s_env, null, new ArrayList<>(),
        OutsideCallParser.parse(s_env, "RYR1\tc.97A>G/[c.97A>G + c.152C>A + c.418G>A]"), null);
    ReportContext reportContext = new ReportContext(s_env, phenotyper, null);

    String json = """
        {
          "guideline": {
            "id": "test:guideline",
            "name": "Test Guideline for testdrug",
            "objCls": "Guideline Annotation",
            "source": "CPIC",
            "recommendation": true,
            "relatedChemicals": [{"id": "test:drug", "name": "testdrug"}],
            "relatedGenes": [{"id": "test:gene", "name": "RYR1", "symbol": "RYR1"}]
          },
          "recommendations": [],
          "citations": [],
          "url": "https://example.com/test-guideline"
        }
        """;
    GuidelinePackage guidelinePackage = DataSerializer.GSON.fromJson(json, GuidelinePackage.class);
    DrugReport drugReport = new DrugReport("testdrug", List.of(guidelinePackage), reportContext);

    Report report = new Report(PrescribingGuidanceSource.CPIC_GUIDELINE, drugReport);
    assertFalse(report.isMatched());
    assertFalse(report.isUnmatchedDpydInferred(), "RYR1 should not trigger the DPYD-specific inferred marker");
    assertTrue(report.isUnmatchedInferred(), "RYR1 should trigger the generic inferred marker");
  }
}
