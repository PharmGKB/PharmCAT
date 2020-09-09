package org.pharmgkb.pharmcat.util;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;


/**
 * This is a JUnit test for {@link CliUtils}.
 *
 * @author Mark Woon
 */
class CliUtilsTest {

  @Test
  void testVersion() throws Exception {

    String version = CliUtils.getVersion();
    assertNotNull(version);
    // we already have tags, so this should never happen if we're running from within a checked-out repo
    assertNotEquals("development", version);
    assertNotEquals("master", version);
  }
}
