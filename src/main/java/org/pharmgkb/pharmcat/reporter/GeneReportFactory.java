package org.pharmgkb.pharmcat.reporter;

import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReportUgt1a1;


/**
 * Factory class for instantiating the right {@link GeneReport} object based on the supplied data
 *
 * @author Ryan Whaley
 */
public class GeneReportFactory {

  /**
   * Make a standard {@link GeneReport} object
   * @param symbol the gene symbol
   * @return an initialized, unpopulated GeneReport
   */
  @Nonnull
  public static GeneReport newReport(@Nonnull String symbol) {
    return new GeneReport(symbol);
  }

  /**
   * Make a specific {@link GeneReport} object based on information found in the GeneCall. Could be a base GeneCall
   * or a gene-specific extended GeneReport class
   * @param geneCall the {@link GeneCall} object to create a report for
   * @return an initialized, unpopulated GeneReport
   */
  @Nonnull
  public static GeneReport newReport(@Nonnull GeneCall geneCall) {
    GeneReport report;

    switch (geneCall.getGene()) {
      case "UGT1A1":
        if (!geneCall.isPhased() || geneCall.getDiplotypes().size() != 1) {
          report = new GeneReportUgt1a1(geneCall.getGene());
        } else {
          report = new GeneReport(geneCall.getGene());
        }
        break;
      default:
        report = new GeneReport(geneCall.getGene());
    }

    return report;
  }
}
