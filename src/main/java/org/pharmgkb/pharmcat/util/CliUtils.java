package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.jar.Attributes;
import java.util.jar.Manifest;
import javax.annotation.Nonnull;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.pharmcat.PharmCAT;


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

  /**
   * Gets the currently tagged version based on the Jar manifest, the current git repo tag, or a generic 
   * "development" version as a fallback when neither of those are available.
   *
   * @return a String of the PharmCAT version
   * @throws IOException can occur from reading the manifest
   */
  public static String getVersion() throws IOException {
    Class clazz = PharmCAT.class;
    String className = clazz.getSimpleName() + ".class";
    String classPath = clazz.getResource(className).toString();
    if (!classPath.startsWith("jar")) {
      // Class not from JAR
      try (StringWriter writer = new StringWriter()) {
        IOUtils.copy(Runtime.getRuntime().exec("git describe --tags").getInputStream(), writer,
            Charset.defaultCharset());
        String gitVersion = writer.toString();
        if (StringUtils.isNotBlank(gitVersion)) {
          return gitVersion;
        }
      } catch (Exception e) {
        System.err.println("Error reading git version");
      }
      return "development";
    }
    String manifestPath = classPath.substring(0, classPath.lastIndexOf("!") + 1) +
        "/META-INF/MANIFEST.MF";
    try (InputStream input = new URL(manifestPath).openStream()) {
      Manifest manifest = new Manifest(input);
      Attributes attr = manifest.getMainAttributes();
      return attr.getValue("Implementation-Version");
    }
  }
}
