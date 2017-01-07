package org.pharmgkb.pharmcat.util;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import javax.annotation.Nonnull;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.io.util.CliHelper;


/**
 * CLI utility methods.
 *
 * @author Mark Woon
 */
public class CliUtils {


  /**
   * Looks for the pharmcat.properties file in the user's home directory.
   *
   * @throws IllegalStateException if properties file cannot be found
   */
  public static @Nonnull Path getPropsFile(@Nonnull CliHelper cliHelper, @Nonnull String propArgKey) {

    if (cliHelper.hasOption(propArgKey)) {
      return cliHelper.getValidFile(propArgKey, true);
    }

    String home = StringUtils.stripToNull(System.getenv("HOME"));
    if (home == null) {
      home = System.getProperty("user.home");
    }
    if (home == null) {
      throw new IllegalStateException("Cannot determine home directory.  Please specify property file.");
    }
    Path propsFile = Paths.get(home).resolve("pharmcat.properties");
    System.out.println("Looking for properties file in: " + propsFile);
    if (!Files.exists(propsFile)) {
      throw new IllegalStateException("Cannot find " + propsFile);
    }
    if (!Files.isRegularFile(propsFile)) {
      throw new IllegalStateException("Not a file: " + propsFile);
    }
    return propsFile;
  }
}
