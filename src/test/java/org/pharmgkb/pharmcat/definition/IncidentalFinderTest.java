package org.pharmgkb.pharmcat.definition;

import org.junit.Test;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;


/**
 * test the IncidentalFinder to make sure it's loading properly and querying exptected results.
 *
 * @author Ryan Whaley
 */
public class IncidentalFinderTest {

  @Test
  public void testConstruct() throws Exception {
    IncidentalFinder incidentalFinder = new IncidentalFinder();

    Haplotype h1 = new Haplotype("CFTR", "R334W");
    assertTrue(incidentalFinder.isFinding(h1));

    Haplotype h2 = new Haplotype("CFTR", "foo");
    assertFalse(incidentalFinder.isFinding(h2));
  }
}
