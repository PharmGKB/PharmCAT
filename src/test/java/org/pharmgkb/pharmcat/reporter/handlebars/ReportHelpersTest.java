package org.pharmgkb.pharmcat.reporter.handlebars;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * This is a JUnit test for {@link ReportHelpers}.
 *
 * @author Mark Woon
 */
class ReportHelpersTest {


  @Test
  void testSanitizeCssSelector() {
    assertEquals("ABC", ReportHelpers.sanitizeCssSelector("__ABC__"));
    assertEquals("ABC", ReportHelpers.sanitizeCssSelector("(ABC)!"));
    assertEquals("A_B-C", ReportHelpers.sanitizeCssSelector("A B - C"));
    assertEquals("MT-RNR1", ReportHelpers.sanitizeCssSelector("MT-RNR1"));
    assertEquals("IFNL3_4", ReportHelpers.sanitizeCssSelector("IFNL3/4"));
  }
}
