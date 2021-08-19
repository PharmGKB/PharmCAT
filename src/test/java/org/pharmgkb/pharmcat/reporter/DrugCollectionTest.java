package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.Set;
import java.util.stream.Collectors;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.assertEquals;


class DrugCollectionTest {

  /**
   * Tests the count of the absolute list of drugs coming in from CPIC. If this assertion fails that means the list of
   * drugs has changed and if you're not expecting that then there may be a problem with the {@link DataManager} or with
   * CPIC.
   */
  @Test
  void testLoad() throws IOException {
    DrugCollection drugCollection = new DrugCollection();
    assertEquals(66, drugCollection.size());
  }

  /**
   * Tests the "reportable" drugs. These are drugs that are NOT associated with any IGNORED_GENES as specified in
   * {@link GeneReport}. If this count changes and you're not expecting it then there may be a porblem with the
   * {@link DataManager}, CPIC datasource, or with how the {@link DrugCollection#listReportable()} method is working.
   */
  @Test
  void testLoadReportable() throws IOException {
    DrugCollection drugCollection = new DrugCollection();
    assertEquals(53, drugCollection.listReportable().size());

    Set<String> reportableDrugNames = drugCollection.listReportable().stream().map((drug) -> drug.getDrugName()).collect(Collectors.toSet());
    Set<String> drugNames = drugCollection.list().stream().map((drug) -> drug.getDrugName()).collect(Collectors.toSet());

    drugNames.removeAll(reportableDrugNames);
    System.out.println("Non-reportable drugs: " + String.join("; ", drugNames));
  }
}
