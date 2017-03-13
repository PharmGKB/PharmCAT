package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Properties;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.SheetsHelper;


/**
 * This class manages allele definition files.
 *
 * @author Mark Woon
 */
public class DefinitionManager {
  public static final Path DEFAULT_DEFINITION_DIR =
      PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles");
  public static final String EXEMPTIONS_FILE_NAME = "exemptions.json";
  public static final String MESSAGES_JSON_FILE = "messages.json";
  private static final Gson sf_gson = new GsonBuilder().serializeNulls().excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting().create();
  private String m_googleUser;
  private String m_googleKey;


  private DefinitionManager(Path propertyFile) throws IOException {

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
          .addOption("m", "messages-dir", "directory to write messages to", false, "m")
          .addOption("d", "download", "download curated allele definition files");

      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      Path propsFile = CliUtils.getPropsFile(cliHelper, "p");
      Path downloadDir = cliHelper.getValidDirectory("in", true);
      Path generatedDir = cliHelper.getValidDirectory("out", true);
      Path messagesDir = null;
      if (cliHelper.hasOption("m")) {
        messagesDir = cliHelper.getValidDirectory("m", true);
      }

      DefinitionManager manager = new DefinitionManager(propsFile);
      if (cliHelper.hasOption("d")) {
        manager.download(downloadDir, generatedDir, messagesDir);
      }
      manager.transform(downloadDir, generatedDir);

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  private void download(@Nonnull Path downloadDir, @Nonnull Path generatedDir, @Nullable Path messagesDir)
      throws Exception {

    SheetsHelper sh = new SheetsHelper(m_googleUser, m_googleKey);
    sh.downloadAlleleDefinitions(downloadDir);
    Path tmpFile = Files.createTempFile("exemptions", ".tsv");
    try {
      sh.downloadAlleleExemptionsFile(tmpFile);
      transformExemptions(tmpFile, generatedDir.resolve(EXEMPTIONS_FILE_NAME));
    } finally {
      Files.deleteIfExists(tmpFile);
    }
    if (messagesDir != null) {
      tmpFile = Files.createTempFile("messages", ".tsv");
      try {
        Path jsonFile = messagesDir.resolve(MESSAGES_JSON_FILE);
        sh.downloadMessagesFile(tmpFile);
        List<MessageAnnotation> messages = Files.lines(tmpFile)
            .skip(1) // skip the header
            .map(MessageAnnotation::new)
            .collect(Collectors.toList());

        try (BufferedWriter writer = Files.newBufferedWriter(jsonFile, StandardCharsets.UTF_8)) {
          writer.write(sf_gson.toJson(messages));
        }
      } finally {
        Files.deleteIfExists(tmpFile);
      }
    }
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private void transform(Path downloadDir, Path outDir) throws IOException {

    GeneratedDefinitionSerializer serializer = new GeneratedDefinitionSerializer();
    try (DirectoryStream<Path> files = Files.newDirectoryStream(downloadDir, f -> f.toString().endsWith("_translation.tsv"))) {
      for (Path file : files) {
        System.out.println("Parsing " + file);
        CuratedDefinitionParser parser = new CuratedDefinitionParser(file);

        DefinitionFile definitionFile = parser.parse();
        if (!parser.getWarnings().isEmpty()) {
          System.out.println("Warnings for " + file);
          parser.getWarnings()
              .forEach(System.out::println);
        }

        Path jsonFile = outDir.resolve(PathUtils.getBaseFilename(file) + ".json");
        serializer.serializeToJson(definitionFile, jsonFile);
        System.out.println("Wrote " + jsonFile);
        System.out.println();
      }
    }
  }

  private void transformExemptions(@Nonnull Path tsvFile, @Nonnull Path jsonFile) throws IOException {
    GeneratedDefinitionSerializer serializer = new GeneratedDefinitionSerializer();
    serializer.serializeToJson(serializer.deserializeExemptionsFromTsv(tsvFile), jsonFile);
  }
}
