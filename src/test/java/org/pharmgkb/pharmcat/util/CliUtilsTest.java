package org.pharmgkb.pharmcat.util;

import org.junit.Test;

import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;


/**
 * This is a JUnit test for {@link CliUtils}.
 *
 * @author Mark Woon
 */
public class CliUtilsTest {

  @Test
  public void testVersion() throws Exception {

    String version = CliUtils.getVersion();
    assertNotNull(version);
    // we already have tags, so this should never happen if we're running from within a checked-out repo
    assertNotEquals("development", version);
  }
}
