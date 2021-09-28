package org.pharmgkb.pharmcat.util;

import java.util.regex.Pattern;
import org.apache.commons.lang3.StringUtils;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;


/**
 * This is a JUnit test for {@link CliUtils}.
 *
 * @author Mark Woon
 */
class CliUtilsTest {
  private static final Pattern sf_versionPattern = Pattern.compile("^v\\d+\\.\\d+\\.\\d+(-\\d+-[a-z0-9]+)?$");

  @Test
  void testVersion() throws Exception {

    String version = StringUtils.stripToNull(CliUtils.getVersion());
    assertNotNull(version);
    System.out.println("Version: [" + version + "]");
    String githubAction = System.getenv("GITHUB_ACTION");
    if (githubAction != null) {
      System.out.println("GITHUB_EVENT_NAME: " + System.getenv("GITHUB_EVENT_NAME"));
      if ("published".equalsIgnoreCase(System.getenv("GITHUB_EVENT_NAME"))) {
        assertTrue(sf_versionPattern.matcher(version).matches());
      } else {
        // when running via GH Actions, won't get anything get via git on push event
        assertEquals("development", version);
      }
    } else {
      assertTrue(sf_versionPattern.matcher(version).matches());
    }
  }
}
