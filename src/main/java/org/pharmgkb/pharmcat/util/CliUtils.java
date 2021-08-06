package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.net.URL;
import java.nio.charset.Charset;
import java.util.Optional;
import java.util.jar.Attributes;
import java.util.jar.Manifest;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.PharmCAT;


/**
 * CLI utility methods.
 *
 * @author Mark Woon
 */
public class CliUtils {


  /**
   * Gets the currently tagged version based on the Jar manifest, the current git repo tag, or a generic 
   * "development" version as a fallback when neither of those are available.
   *
   * @return a String of the PharmCAT version
   * @throws IOException can occur from reading the manifest
   */
  public static String getVersion() throws IOException {
    Class<PharmCAT> clazz = PharmCAT.class;
    String className = clazz.getSimpleName() + ".class";
    String classPath = Optional.ofNullable(clazz.getResource(className))
        .orElseThrow(() -> new IOException("Unable to find class " + className))
        .toString();
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
