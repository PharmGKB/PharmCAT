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
  private static Path s_testOutputDir = getDefaultTestOutputDir();
  private static boolean s_saveTestOutput = "true".equals(System.getenv("PHARMCAT_SAVE_TEST_OUTPUT"));


  private TestUtils() {
  }

  public static void setTestOutputDir(Path outputDir) {
    s_testOutputDir = outputDir;
  }

  public static Path getTestOutputDir() {
    return s_testOutputDir;
  }


  public static boolean isSaveTestOutput() {
    return s_saveTestOutput;
  }

  public static void setSaveTestOutput(boolean saveTestOutput) {
    s_saveTestOutput = saveTestOutput;
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



  private static Path getDefaultTestOutputDir() {
    Path outputDir;
    if (System.getProperty("PHARMCAT_TEST_DIR") != null) {
      outputDir = Paths.get(System.getProperty("PHARMCAT_TEST_DIR"));
    } else {
      outputDir = Paths.get("out");
      if (!outputDir.isAbsolute()) {
        outputDir = Paths.get(System.getProperty("user.dir")).resolve("tmp/test_output");
      } else {
        outputDir = outputDir.resolve("test_output");
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


  /**
   * Gets the output directory for the given test.
   * Directory is guaranteed to exist.
   *
   * @param deleteIfExist if directory exists, it will be deleted and re-created
   */
  public static Path getTestOutputDir(TestInfo testInfo, boolean deleteIfExist) throws IOException {
    Path classOutputDir;
    if (testInfo.getTestClass().isPresent()) {
      classOutputDir = s_testOutputDir.resolve(testInfo.getTestClass().get().getSimpleName());
    } else {
      classOutputDir = s_testOutputDir;
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
    if (!Files.exists(dir)) {
      Files.createDirectories(dir);
    }
    return dir;
  }

  public static Path getTestOutputDir(Class testClass, boolean deleteIfExist) throws IOException {
    Path dir = s_testOutputDir.resolve(testClass.getSimpleName());
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


  /**
   * Creates a temporary test file based on {@code testInfo}, with the specified suffix.
   * This file name will be change from one test to the next.
   */
  public static Path createTempFile(TestInfo testInfo, String suffix) throws IOException {
    return createTempFile(getTestOutputDir(testInfo, false), getTestName(testInfo), suffix);
  }

  private static Path createTempFile(Path dir, String prefix, String suffix) throws IOException {
    if (!Files.isDirectory(dir)) {
      Files.createDirectories(dir);
    }
    Path file = Files.createTempFile(dir, prefix, suffix);
    if (!s_saveTestOutput) {
      file.toFile().deleteOnExit();
    }
    return file;
  }

  public static Path createTempDirectory(String prefix) throws IOException {
    return Files.createTempDirectory(s_testOutputDir, prefix);
  }


  /**
   * Creates a test file based on {@code testInfo}, with the specified suffix.
   * This file name will be constant from one test to the next.
   */
  public static Path createTestFile(TestInfo testInfo, String suffix) throws IOException {
    return createTestFile(getTestOutputDir(testInfo, false), getTestName(testInfo) + suffix);
  }

  public static Path createTestFile(Class testClass, String filename) throws IOException {
    return createTestFile(getTestOutputDir(testClass, false), filename);
  }

  /**
   * Creates a test file with the given filename.
   * This file name will be constant from one test to the next.
   */
  private static Path createTestFile(Path dir, String filename) throws IOException {
    if (!Files.isDirectory(dir)) {
      Files.createDirectories(dir);
    }
    Path file = dir.resolve(filename);
    if (!s_saveTestOutput) {
      file.toFile().deleteOnExit();
    }
    return file;
  }

  public static void deleteTestFiles(Path... files) throws IOException {
    if (!isSaveTestOutput()) {
      for (Path file : files) {
        Files.deleteIfExists(file);
      }
    }
  }
}
