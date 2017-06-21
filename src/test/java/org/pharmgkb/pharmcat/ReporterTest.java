package org.pharmgkb.pharmcat;

import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.Test;
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
  public void testCypc2c9VariantPassthrough() {
    String vcfData = VcfTestUtils.writeVcf(new String[]{
        "CYP2C9/s1s1.vcf"
    });

    try {
      Path tempOutDir = Files.createTempDirectory("ReporterTest");
      Path sampleVcf = Files.createTempFile(tempOutDir, "sample",".vcf");
      Files.write(sampleVcf, vcfData.getBytes());

      PharmCAT pharmcat = new PharmCAT(tempOutDir, null, null);
      pharmcat.execute(sampleVcf, null, null);

      Reporter reporter = pharmcat.getReporter();
      GeneReport geneReport = reporter.getContext().getGeneReport("CYP2C9");

      assertNotNull(geneReport);
      assertNotNull(geneReport.getVariantReports());


      assertTrue(
          "Exemption variant not included in gene report",
          geneReport.getVariantReports().stream()
              .anyMatch(r -> r.getDbSnpId() != null && r.getDbSnpId().equals("rs12777823"))
      );
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
