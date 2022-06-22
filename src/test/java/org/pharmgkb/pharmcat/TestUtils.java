package org.pharmgkb.pharmcat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import org.apache.commons.io.FileUtils;
import org.junit.jupiter.api.TestInfo;


/**
 * Helper functions for working with tests.
 *
 * @author Mark Woon
 */
public class TestUtils {
  public static final Path TEST_OUTPUT_DIR = getTestOutputDir();

  private TestUtils() {
  }


  public static String getTestName(TestInfo testInfo) {
    return testInfo.getDisplayName().replace("(TestInfo)", "");
  }


  /**
   * Checks if test is running in a continuous integration environment.
   * This is determined based on the `CI` envinorment variable on GH Actions:
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
    Path outputDir;
    if (System.getProperty("PHARMCAT_TEST_DIR") != null) {
      outputDir = Paths.get(System.getProperty("PHARMCAT_TEST_DIR"));
    } else {
      outputDir = Paths.get("out");
      if (!outputDir.isAbsolute()) {
        outputDir = Paths.get(System.getProperty("user.dir")).resolve("out");
      }
    }
    if (!Files.isDirectory(outputDir)) {
      try {
        Files.createDirectories(outputDir);
      } catch (IOException ex) {
        throw new RuntimeException("Error creating test directory: " + outputDir);
      }
    }
    return outputDir;
  }


  public static Path getTestOutputDir(TestInfo testInfo, boolean deleteIfExist) throws IOException {
    Path classOutputDir;
    if (testInfo.getTestClass().isPresent()) {
      classOutputDir = TEST_OUTPUT_DIR.resolve(testInfo.getTestClass().get().getSimpleName());
    } else {
      classOutputDir = TEST_OUTPUT_DIR;
    }
    Path dir = classOutputDir.resolve(getTestName(testInfo));
    if (Files.exists(dir)) {
      if (Files.isDirectory(dir)) {
        if (deleteIfExist) {
          FileUtils.deleteDirectory(dir.toFile());
        }
      } else {
        throw new RuntimeException("Not a directory: " + dir);
      }
    }
    return dir;
  }

  public static Path createTempFile(String prefix, String suffix) throws IOException {
    return createTempFile(TEST_OUTPUT_DIR, prefix, suffix);
  }

  public static Path createTempFile(TestInfo testInfo, String suffix) throws IOException {
    return createTempFile(getTestOutputDir(testInfo, false), getTestName(testInfo), suffix);
  }

  private static Path createTempFile(Path dir, String prefix, String suffix) throws IOException {
    if (!Files.isDirectory(dir)) {
      Files.createDirectories(dir);
    }
    Path file = Files.createTempFile(dir, prefix, suffix);
    file.toFile().deleteOnExit();
    return file;
  }

  public static Path createTempDirectory(String prefix) throws IOException {
    return Files.createTempDirectory(TEST_OUTPUT_DIR, prefix);
  }
}
