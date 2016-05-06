package org.pharmgkb.pharmcat.definition;

import java.io.File;
import java.nio.file.Path;
import org.junit.Test;
import org.pharmgkb.pharmcat.TestUtil;

import static org.junit.Assert.*;


/**
 * @author Ryan Whaley
 */
public class GuidelineFileSetTest {

  @Test
  public void test() throws Exception {
    Path testGuidelineDir = TestUtil.getFile("org/pharmgkb/pharmcat/annotations");

    assertNotNull(testGuidelineDir);
    assertTrue(testGuidelineDir.toFile().exists());
    assertTrue(testGuidelineDir.toFile().isDirectory());

    File[] files = testGuidelineDir.toFile().listFiles();
    assertNotNull(files);
    assertEquals(30, files.length);


  }
}
