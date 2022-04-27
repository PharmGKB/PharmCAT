package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.util.Set;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.DataManager;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static org.junit.jupiter.api.Assertions.assertEquals;


class DrugCollectionTest {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  /**
   * Tests the count of the absolute list of drugs coming in from CPIC. If this assertion fails that means the list of
   * drugs has changed and if you're not expecting that then there may be a problem with the {@link DataManager} or with
   * CPIC.
   */
  @Test
  void testLoad() throws IOException {
    DrugCollection drugCollection = new DrugCollection();
    assertEquals(72, drugCollection.list().stream().filter(DrugCollection.ONLY_CPIC).count());
    assertEquals(55, drugCollection.list().stream().filter(DrugCollection.ONLY_DPWG).count());
  }

  /**
   * Tests the "reportable" drugs. These are drugs that are NOT associated with any IGNORED_GENES as specified in
   * {@link GeneReport}. If this count changes and you're not expecting it then there may be a porblem with the
   * {@link DataManager}, CPIC datasource, or with how the {@link DrugCollection#listReportable()} method is working.
   */
  @Test
  void testLoadReportable() throws IOException {
    DrugCollection drugCollection = new DrugCollection();
    assertEquals(127, drugCollection.listReportable().size());

    Set<String> reportableDrugNames = drugCollection.listReportable().stream().map(Drug::getDrugName).collect(Collectors.toSet());
    Set<String> drugNames = drugCollection.list().stream().map(Drug::getDrugName).collect(Collectors.toSet());

    drugNames.removeAll(reportableDrugNames);
    sf_logger.info("Non-reportable drugs: " + String.join("; ", drugNames));
  }
}
