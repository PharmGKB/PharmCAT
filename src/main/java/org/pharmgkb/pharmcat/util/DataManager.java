package org.pharmgkb.pharmcat.util;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;
import java.lang.invoke.MethodHandles;
import java.net.URL;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableSet;
import com.google.gson.Gson;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.pharmgkb.common.io.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.reporter.DrugCollection;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;


/**
 * This class manages external resources (e.g. allele definition files, dosing guideline annotations)
 *
 * @author Mark Woon
 */
public class DataManager {
  public static final Path DEFAULT_DEFINITION_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles");
  public static final String EXEMPTIONS_JSON_FILE_NAME = "exemptions.json";
  private static final String MESSAGES_JSON_FILE_NAME = "messages.json";
  private static final String SUMMARY_REPORT = "summary.md";
  private static final Set<String> PREFER_OUTSIDE_CALL = ImmutableSet.of("CYP2D6");
  private final DataSerializer m_dataSerializer = new DataSerializer();
  private final boolean m_verbose;


  private DataManager(boolean verbose) {
    m_verbose = verbose;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("dl", "download-dir", "directory to save downloaded files", false, "dl")
          .addOption("a", "alleles-dir", "directory to save generated allele definition files", true, "a")
          .addOption("m", "messages-dir", "directory to write messages to", true, "m")
          .addOption("g", "drugs-dir", "directory to save drug data to", false, "g")
          .addOption("d", "documenation-dir", "directory to save documentation to", false, "documentation-path")
          .addOption("sd", "skip-download", "skip downloading")
          .addOption("sa", "skip-alleles", "skip alleles")
          .addOption("sm", "skip-messages", "skip messages")
          .addOption("sg", "skip-drugs", "skip drugs");


      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      boolean skipDownload = cliHelper.hasOption("sd");
      boolean skipDrugs = cliHelper.hasOption("sg");
      boolean skipAlleles = cliHelper.hasOption("sa");

      Path downloadDir;
      if (cliHelper.hasOption("dl")) {
        downloadDir = cliHelper.getValidDirectory("dl", true);
      } else {
        downloadDir = Files.createTempDirectory("pharmcat");
        downloadDir.toFile().deleteOnExit();
      }

      try {
        Path allelesDir = cliHelper.getValidDirectory("a", true);
        Path messageDir = cliHelper.getValidDirectory("m", true);

        Path drugsDir = null;
        if (cliHelper.hasOption("g")) {
          drugsDir = cliHelper.getValidDirectory("g", true);
        }

        Path exemptionsTsv = downloadDir.resolve("exemptions.tsv");
        Path exemptionsJson = allelesDir.resolve(EXEMPTIONS_JSON_FILE_NAME);
        Path messagesTsv = downloadDir.resolve("messages.tsv");
        Path messagesJson = messageDir.resolve(MESSAGES_JSON_FILE_NAME);

        DataManager manager = new DataManager(cliHelper.isVerbose());
        if (!skipDownload) {
          manager.download(downloadDir, exemptionsTsv, messagesTsv);
        }

        DrugCollection drugs = manager.saveDrugData(drugsDir, skipDrugs);

        Set<String> genes;
        if (!skipAlleles) {
          Map<String, DefinitionExemption> exemptionsMap = manager.transformExemptions(exemptionsTsv, exemptionsJson);
          // if we're loading new gene data, use the list of genes from the new data
          genes = manager.transformAlleleDefinitions(downloadDir, allelesDir, exemptionsMap);
        } else {
          // if we're sipping new gene data, then use the existing list of genes
          genes = new DefinitionReader().getGenes();
        }

        if (!cliHelper.hasOption("sm")) {
          manager.transformMessages(messagesTsv, messagesJson);
        }

        if (cliHelper.hasOption("d")) {
          Path docsDir = cliHelper.getValidDirectory("d", true);
          manager.writeSummary(docsDir, genes, drugs);
        }

      } finally {
        if (!cliHelper.hasOption("dl")) {
          FileUtils.deleteQuietly(downloadDir.toFile());
        }
      }

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  private void download(Path downloadDir, Path exemptionsFile, Path messagesFile) throws Exception {

    // download allele definitions
    FileUtils.copyURLToFile(
        new URL("https://files.cpicpgx.org/data/report/current/allele_definitions.json"),
        downloadDir.resolve("allele_definitions.json").toFile());

    String urlFmt = "https://docs.google.com/spreadsheets/d/%s/export?format=tsv";
    // download exemptions and messages
    FileUtils.copyURLToFile(new URL(String.format(urlFmt, "1xHvvXQIMv3xbqNhuN7zG6WP4DB7lpQDmLvz18w-u_lk")),
        exemptionsFile.toFile());
    FileUtils.copyURLToFile(new URL(String.format(urlFmt, "1MkWV6TlJTnw-KRNWeylyUJAUocCgupcJLlmV2fRdtcM")),
        messagesFile.toFile());
  }

  private DrugCollection saveDrugData(Path drugsDir, boolean skip) throws Exception {
    DrugCollection existingDrugs = new DrugCollection();
    if (drugsDir == null || skip) return existingDrugs;

    Path drugsPath = drugsDir.resolve("drugs.json");
    FileUtils.copyURLToFile(
        new URL("https://files.cpicpgx.org/data/report/current/drugs.json"),
        drugsPath.toFile());
    System.out.println();
    System.out.println("Saving messages to " + drugsPath);

    Set<String> existingDrugNames = existingDrugs.list().stream()
        .map(Drug::getDrugName).collect(Collectors.toSet());
    try (InputStream is = Files.newInputStream(drugsPath)) {
      DrugCollection newDrugs = new DrugCollection(is);
      Set<String> newDrugNames = newDrugs.list().stream()
          .map(Drug::getDrugName).collect(Collectors.toSet());

      Set<String> removedDrugs = new HashSet<>(existingDrugNames);
      removedDrugs.removeAll(newDrugNames);
      if (removedDrugs.size() > 0) {
        System.out.println("Removed drugs: " + String.join("; ", removedDrugs));
      }

      Set<String> addedDrugs = new HashSet<>(newDrugNames);
      addedDrugs.removeAll(existingDrugNames);
      if (addedDrugs.size() > 0) {
        System.out.println("New drugs: " + String.join("; ", addedDrugs));
      }
      return newDrugs;
    }
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private Set<String> transformAlleleDefinitions(Path downloadDir, Path definitionsDir,
      Map<String, DefinitionExemption> exemptionsMap) throws Exception {

    Path definitionsFile = downloadDir.resolve("allele_definitions.json");
    if (!Files.exists(definitionsFile)) {
      throw new IOException("Cannot find alleles definitions (" + definitionsFile + ")");
    }
    String json = FileUtils.readFileToString(definitionsFile.toFile(), Charsets.UTF_8);
    Gson gson = new Gson();
    Map<String, DefinitionFile> definitionFileMap = new HashMap<>();
    for (DefinitionFile df : gson.fromJson(json, DefinitionFile[].class)) {
      // TODO: remove after v1 is released
      if (df.getGeneSymbol().equals("MT-RNR1")) continue;

      definitionFileMap.put(df.getGeneSymbol(), df);
    }

    System.out.println();
    System.out.println("Saving allele definitions in " + definitionsDir.toString());
    Set<String> currentFiles = Files.list(definitionsDir)
        .map(PathUtils::getFilename)
        .filter(f -> f.endsWith("_translation.json"))
        .collect(Collectors.toSet());

    // transform data as necessary
    for (String gene : exemptionsMap.keySet()) {
      DefinitionExemption exemption = exemptionsMap.get(gene);
      if (exemption.getIgnoredPositions().size() > 0) {
        System.out.println("Removing ignored positions in " + gene + "...");
        DefinitionFile definitionFile = definitionFileMap.get(gene);
        if (definitionFile == null) {
          throw new Exception("No definition for " + gene);
        }
        definitionFile.removeIgnoredPositions(exemption);
      }
    }
    fixCyp2c19(definitionFileMap.get("CYP2C19"));

    // output file
    for (String gene : definitionFileMap.keySet()) {
      DefinitionFile definitionFile = definitionFileMap.get(gene);
      Path jsonFile = definitionsDir.resolve(gene + "_translation.json");
      m_dataSerializer.serializeToJson(definitionFile, jsonFile);
      if (m_verbose) {
        System.out.println("Wrote " + jsonFile);
      }
      if (!currentFiles.remove(gene + "_translation.json")) {
        System.out.println("New gene: " + gene);
      }
    }

    deleteObsoleteFiles(definitionsDir, currentFiles);
    return definitionFileMap.keySet();
  }

  private void fixCyp2c19(DefinitionFile definitionFile) {
    Preconditions.checkNotNull(definitionFile);
    NamedAllele star1 = definitionFile.getNamedAlleles().stream()
        .filter(na -> na.getName().equals("*1"))
        .findFirst()
        .orElseThrow(() -> new IllegalStateException("Cannot find CYP2C19*1"));
    NamedAllele star38 = definitionFile.getNamedAlleles().stream()
        .filter(na -> na.getName().equals("*38"))
        .findFirst()
        .orElseThrow(() -> new IllegalStateException("Cannot find CYP2C19*38"));
    for (int x = 0; x < star38.getAlleles().length; x += 1) {
      if (star1.getAlleles()[x] == null) {
        star1.getAlleles()[x] = star38.getAlleles()[x];
      }
    }
  }


  private Map<String, DefinitionExemption> transformExemptions(Path tsvFile, Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving exemptions to " + jsonFile.toString());
    Set<DefinitionExemption> exemptions = m_dataSerializer.deserializeExemptionsFromTsv(tsvFile);
    m_dataSerializer.serializeToJson(exemptions, jsonFile);

    Map<String, DefinitionExemption> exemptionsMap = new HashMap<>();
    exemptions.forEach(exemption -> exemptionsMap.put(exemption.getGene(), exemption));
    return exemptionsMap;
  }

  private void transformMessages(Path tsvFile, Path jsonFile) throws IOException {

    System.out.println();
    System.out.println("Saving messages to " + jsonFile.toString());
    m_dataSerializer.serializeToJson(m_dataSerializer.deserializeMessagesFromTsv(tsvFile), jsonFile);
  }


  private void deleteObsoleteFiles(Path dir, Set<String> obsoleteFilenames) {

    for (String filename : obsoleteFilenames) {
      Path file = dir.resolve(filename);
      System.out.println("Deleting obsolete file: " + file);
      FileUtils.deleteQuietly(file.toFile());
    }
  }

  private void writeSummary(Path documenationDir, Set<String> genes, DrugCollection drugs) throws IOException {
    try (
        InputStream mdStream = getClass().getResourceAsStream("summary.md");
        StringWriter templateWriter = new StringWriter()
    ) {
      if (mdStream == null) {
        throw new RuntimeException("No summary report templater found");
      }
      IOUtils.copy(mdStream, templateWriter, Charset.defaultCharset());
      String mdTemplate = templateWriter.toString();

      File summaryFile = documenationDir.resolve(SUMMARY_REPORT).toFile();
      try (FileWriter fw = new FileWriter(summaryFile)) {
        String matcherGeneList = genes.stream()
            .filter(g -> !PREFER_OUTSIDE_CALL.contains(g))
            .sorted().map(g -> "- " + g)
            .collect(Collectors.joining("\n"));
        String outsideGeneList = genes.stream()
            .filter(PREFER_OUTSIDE_CALL::contains)
            .sorted().map(g -> "- " + g)
            .collect(Collectors.joining("\n"));
        String drugList = drugs.listReportable().stream().map(d -> "- " + d.getDrugName()).collect(Collectors.joining("\n"));
        IOUtils.write(String.format(mdTemplate, matcherGeneList, outsideGeneList, drugList), fw);
      }
      System.out.println("\nSaving summary file to " + summaryFile);
    }
  }
}
