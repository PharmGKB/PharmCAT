package org.pharmgkb.pharmcat;

/**
 * Helper functions for working with tests.
 *
 * @author Mark Woon
 */
public class TestUtils {


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
}
