package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Path;
import java.util.Optional;
import java.util.jar.Attributes;
import java.util.jar.Manifest;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.PharmCAT;


/**
 * CLI utility methods.
 *
 * @author Mark Woon
 */
public class CliUtils {
  private static String s_version;


  /**
   * Private constructor for utility class.
   */
  private CliUtils() {
  }


  /**
   * Gets the currently tagged version based on the Jar manifest, the current git repo tag, or a generic 
   * "development" version as a fallback when neither of those are available.
   *
   * @return a String of the PharmCAT version
   * @throws IOException can occur from reading the manifest
   */
  public synchronized static String getVersion() throws IOException {
    if (s_version != null) {
      return s_version;
    }
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
        String gitVersion = StringUtils.strip(writer.toString());
        if (StringUtils.isNotBlank(gitVersion)) {
          s_version = gitVersion;
          return s_version;
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
      s_version = attr.getValue("Implementation-Version");
      return s_version;
    }
  }


  public static Path getOutputDir(CliHelper cliHelper, Path inputFile) throws IOException {
    if (cliHelper.hasOption("o")) {
      return cliHelper.getValidDirectory("o", true);
    } else {
      return inputFile.getParent();
    }
  }

  public static Path getOutputFile(CliHelper cliHelper, Path inputFile, @Nullable String preferredFileArg,
      String defaultSuffix) throws IOException {

    Path preferredFile = null;
    if (preferredFileArg != null && cliHelper.getValue(preferredFileArg) != null) {
      preferredFile = cliHelper.getPath(preferredFileArg);
    }
    return getOutputFile(getOutputDir(cliHelper, inputFile), inputFile, preferredFile, defaultSuffix);
  }


  public static Path getOutputFile(Path outputDir, Path inputFile, @Nullable Path preferredFile,
      String defaultSuffix) {

    if (preferredFile != null) {
      if (preferredFile.isAbsolute()) {
        return preferredFile;
      }
      return outputDir.resolve(preferredFile);
    }

    String baseFilename = FilenameUtils.getBaseName(inputFile.getFileName().toString());
    return outputDir.resolve(baseFilename + defaultSuffix);
  }
}
