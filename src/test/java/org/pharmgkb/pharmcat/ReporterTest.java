package org.pharmgkb.pharmcat;

import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.Test;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.Reporter;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;


/**
 * Test run the Reporter and check generated data
 *
 * @author Ryan Whaley
 */
public class ReporterTest {

  @Test
  public void testCypc2c9VariantPassthrough() throws Exception {

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c9/s1s1.vcf");

    Path tempOutDir = Files.createTempDirectory("ReporterTest");

    PharmCAT pharmcat = new PharmCAT(tempOutDir, null, null);
    pharmcat.execute(vcfFile, null, null);

    Reporter reporter = pharmcat.getReporter();
    GeneReport geneReport = reporter.getContext().getGeneReport("CYP2C9");

    assertNotNull(geneReport);
    assertNotNull(geneReport.getVariantReports());


    assertTrue(
        "Exemption variant not included in gene report",
        geneReport.getVariantReports().stream()
            .anyMatch(r -> r.getDbSnpId() != null && r.getDbSnpId().equals("rs12777823"))
    );
  }
}
