package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Properties;
import com.google.common.base.Preconditions;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.io.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;


/**
 * This class manages allele definitions.
 *
 * @author Mark Woon
 */
public class DefinitionManager {
  private String m_googleUser;
  private String m_googleKey;


  public DefinitionManager(Path propertyFile) throws IOException {

    Properties properties = new Properties();
    try (BufferedReader reader = Files.newBufferedReader(propertyFile)) {
      properties.load(reader);
    }
    m_googleUser = StringUtils.stripToNull((String)properties.get("google.user"));
    Preconditions.checkState(m_googleUser != null, "Missing property: 'google.user");
    m_googleKey = StringUtils.stripToNull((String)properties.get("google.key"));
    Preconditions.checkState(m_googleKey != null, "Missing property: 'google.key");
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("p", "properties-file", "PharmCAT properties file", false, "p")
          .addOption("in", "download-dir", "directory of save curated allele definition files", true, "in")
          .addOption("out", "generated-dir", "directory of save generated allele definition files", true, "out")
          .addOption("d", "download", "download curated allele definition files");

      if (cliHelper.parse(args)) {
        System.exit(1);
      }

      Path propsFile;
      if (cliHelper.hasOption("p")) {
        propsFile = cliHelper.getValidFile("p", true);
      } else {
        String home = StringUtils.stripToNull(System.getenv("HOME"));
        if (home == null) {
          home = System.getProperty("user.home");
        }
        if (home == null) {
          System.out.println("Cannot determine home directory.  Please specify property file.");
          System.exit(1);
        }
        propsFile = Paths.get(home).resolve("pharmcat.properties");
        System.out.println("Looking for properties file in: " + propsFile);
        if (!Files.exists(propsFile)) {
          System.out.println("Cannot find " + propsFile);
          System.exit(1);
        }
        if (!Files.isRegularFile(propsFile)) {
          System.out.println("Not a file: " + propsFile);
          System.exit(1);
        }
      }

      Path downloadDir = cliHelper.getValidDirectory("in", true);
      Path generatedDir = cliHelper.getValidDirectory("out", true);

      DefinitionManager manager = new DefinitionManager(propsFile);
      if (cliHelper.hasOption("d")) {
        manager.download(downloadDir);
      }
      manager.transform(downloadDir, generatedDir);

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  private void download(Path downloadDir) throws Exception {

    SheetsHelper sh = new SheetsHelper(m_googleUser, m_googleKey);
    sh.downloadAlleleDefinitions(downloadDir);
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private void transform(Path downloadDir, Path outDir) throws IOException {

    GeneratedDefinitionSerializer serializer = new GeneratedDefinitionSerializer();
    try (DirectoryStream<Path> files = Files.newDirectoryStream(downloadDir, f -> f.toString().endsWith(".tsv"))) {
      for (Path file : files) {
        System.out.println("Parsing " + file);
        CuratedDefinitionParser parser = new CuratedDefinitionParser(file);

        DefinitionFile definitionFile = parser.parse();
        if (parser.getWarnings().isEmpty()) {
          System.out.println("Warnings for " + file);
          parser.getWarnings().stream()
              .forEach(System.out::println);
        }

        serializer.serializeToTsv(definitionFile, outDir.resolve(definitionFile.getGeneSymbol() + ".tsv"));
        Path jsonFile = outDir.resolve(PathUtils.getBaseFilename(file) + ".json");
        serializer.serializeToJson(definitionFile, jsonFile);
        System.out.println("Wrote " + jsonFile);
        System.out.println();
      }
    }
  }
}
