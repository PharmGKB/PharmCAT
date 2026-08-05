package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.pharmcat.TestUtils;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;


/**
 * JUnit test for {@link ZipUtils}.
 *
 * @author Mark Woon
 */
class ZipUtilsTest {

  @AfterEach
  void deleteDirectory(TestInfo testInfo) throws IOException {
    TestUtils.deleteTestOutputDirectory(testInfo);
    // in case the Zip Slip protection ever regresses, make sure the escaped file doesn't linger between runs
    Files.deleteIfExists(TestUtils.getTestOutputDir(testInfo, false).resolveSibling("evil.txt"));
  }


  @Test
  void unzipExtractsNormalEntries(TestInfo testInfo) throws Exception {
    Path destDir = TestUtils.getTestOutputDir(testInfo, true);
    Path zipFile = destDir.resolve("normal.zip");
    writeZip(zipFile, "subdir/file.txt", "hello");

    ZipUtils.unzip(zipFile, destDir);

    Path extracted = destDir.resolve("subdir").resolve("file.txt");
    assertTrue(Files.isRegularFile(extracted));
    assertEquals("hello", Files.readString(extracted, StandardCharsets.UTF_8));
  }


  /**
   * Regression test for Zip Slip: a zip entry using "../" to escape the destination directory must be rejected,
   * not silently extracted outside of it.
   */
  @Test
  void unzipRejectsZipSlip(TestInfo testInfo) throws Exception {
    Path destDir = TestUtils.getTestOutputDir(testInfo, true);
    Path zipFile = destDir.resolve("evil.zip");
    writeZip(zipFile, "../evil.txt", "pwned");

    assertThrows(IOException.class, () -> ZipUtils.unzip(zipFile, destDir));
    assertFalse(Files.exists(destDir.resolveSibling("evil.txt")));
  }


  private void writeZip(Path zipFile, String entryName, String content) throws IOException {
    try (ZipOutputStream zos = new ZipOutputStream(Files.newOutputStream(zipFile))) {
      zos.putNextEntry(new ZipEntry(entryName));
      zos.write(content.getBytes(StandardCharsets.UTF_8));
      zos.closeEntry();
    }
  }
}
