package org.pharmgkb.pharmcat.reporter;

import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


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
    return new GeneReport(geneCall.getGene());
  }
}
