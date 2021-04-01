package org.pharmgkb.pharmcat.util;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;


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
    String githubAction = System.getenv("GITHUB_ACTION");
    if (githubAction != null) {
      // running via GH Actions, won't get anything via git
      assertEquals("development", version);
    } else {
      // we already have tags, so this should never happen if we're running from within a checked-out repo
      assertNotEquals("development", version);
      assertNotEquals("main", version);
    }
  }
}
