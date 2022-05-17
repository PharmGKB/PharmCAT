package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;


/**
 * Helper functions for working with tests.
 *
 * @author Mark Woon
 */
public class TestUtils {
  public static final Path TEST_OUTPUT_DIR = getTestOutputDir();

  private TestUtils() {
  }


  /**
   * Checks if test is running in a continuous integration environment.
   * This is determiend based on the `CI` envinorment variable on GH Actions:
   * https://docs.github.com/en/actions/learn-github-actions/environment-variables#default-environment-variables
   */
  public static boolean isContinuousIntegration() {
    return Boolean.parseBoolean(System.getenv("CI"));
  }


  public static boolean isIgnorableExternalServiceException(Throwable ex) {
    if (ex.getMessage() != null &&
        ex.getMessage().contains("500 error") &&
        ex.getMessage().contains("Problem with 3rd party API, please retry")) {
      System.out.println("IGNORING: " + ex.getMessage());
      return true;
    }
    return false;
  }


  private static Path getTestOutputDir() {
    Path outputDir = Paths.get("out");
    if (!outputDir.isAbsolute()) {
      outputDir = Paths.get(System.getProperty("user.dir")).resolve("out");
    }
    if (!Files.isDirectory(outputDir)) {
      try {
        Files.createDirectories(outputDir);
      } catch (IOException ex) {
        throw new RuntimeException("Error creating directory: " + outputDir);
      }
    }
    return outputDir;
  }

  public static Path createTempFile(String prefix, String suffix) throws IOException {
    Path file = Files.createTempFile(TEST_OUTPUT_DIR, prefix, suffix);
    file.toFile().deleteOnExit();
    return file;
  }

  public static Path createTempDirectory(String prefix) throws IOException {
    Path dir = Files.createTempDirectory(TEST_OUTPUT_DIR, prefix);
    return dir;
  }
}
