package org.pharmgkb.pharmcat.definition;

import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * test the IncidentalFinder to make sure it's loading properly and querying exptected results.
 *
 * @author Ryan Whaley
 */
class IncidentalFinderTest {

  @Test
  void testConstruct() throws Exception {
    IncidentalFinder incidentalFinder = new IncidentalFinder();

    Haplotype h1 = new Haplotype("CFTR", "R334W");
    assertTrue(incidentalFinder.isFinding(h1));

    Haplotype h2 = new Haplotype("CFTR", "foo");
    assertFalse(incidentalFinder.isFinding(h2));
  }
}
