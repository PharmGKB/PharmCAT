package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.lang.invoke.MethodHandles;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import org.apache.commons.io.FileUtils;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.InternalWrapper;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.DiplotypeRecord;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.model.DataSource;


/**
 * This class manages external resources (e.g. allele definition files, dosing guideline annotations).
 *
 * @author Mark Woon
 */
public class DataManager {
  public static final Path DEFAULT_DEFINITION_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles");
  public static final String EXEMPTIONS_JSON_FILE_NAME = "exemptions.json";
  public static final Path DEFAULT_EXEMPTIONS_FILE = DEFAULT_DEFINITION_DIR.resolve(EXEMPTIONS_JSON_FILE_NAME);
  public static final String POSITIONS_VCF = "pharmcat_positions.vcf";
  public static final String UNIALLELIC_POSITIONS_VCF = "pharmcat_positions.uniallelic.vcf";
  private static final String ALLELES_FILE_NAME = "allele_translations.json";
  private static final String sf_zipFileName = "pharmcat.zip";
  private static final String sf_googleDocUrlFmt = "https://docs.google.com/spreadsheets/d/%s/export?format=tsv";

  private final DataSerializer m_dataSerializer = new DataSerializer();
  private final boolean m_verbose;


  private DataManager(boolean verbose) {
    m_verbose = verbose;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("dl", "download-dir", "directory to save downloaded files", false, "dl")
          .addOption("sdl", "skip-download", "skip download")
          .addOption("a", "alleles-dir", "directory to save generated allele definition files", false, "a")
          .addOption("sa", "skip-alleles", "skip alleles")
          .addOption("g", "guidelines-dir", "directory to save guidelines data to", false, "dir")
          .addOption("sg", "skip-drugs", "skip drugs")
          .addOption("m", "messages-dir", "directory to write messages to", false, "m")
          .addOption("sm", "skip-messages", "skip messages")
          .addOption("p", "phenotypes-dir", "directory to save phenotypes to", false, "p")
          .addOption("sp", "skip-phenotypes", "skip phenotype files")
          .addOption("doc", "documentation-dir", "directory to save documentation to", false, "dir")
          ;


      if (!cliHelper.parse(args)) {
        System.exit(1);
      }

      boolean skipDownload = cliHelper.hasOption("sdl");
      Path downloadDir;
      if (cliHelper.hasOption("dl")) {
        downloadDir = cliHelper.getValidDirectory("dl", true);
        System.out.println("Downloading to " + downloadDir);
      } else {
        if (skipDownload) {
          System.out.println("Cannot skip download without providing download directory containing necessary files");
          System.exit(1);
        }
        downloadDir = Files.createTempDirectory("pharmcat");
        downloadDir.toFile().deleteOnExit();
      }

      try {
        DataManager manager = new DataManager(cliHelper.isVerbose());

        if (!cliHelper.hasOption("sm")) {
          Path messagesTsv = downloadDir.resolve("messages.tsv");
          if (!skipDownload) {
            // download messages
            FileUtils.copyURLToFile(new URL(String.format(sf_googleDocUrlFmt, "1MkWV6TlJTnw-KRNWeylyUJAUocCgupcJLlmV2fRdtcM")),
                messagesTsv.toFile());
          }

          Path messageDir = cliHelper.getValidDirectory("m", true);
          Path messagesJson = messageDir.resolve(MessageHelper.MESSAGES_JSON_FILE_NAME);
          manager.transformMessages(messagesTsv, messagesJson);
        }

        boolean skipGuidelines = cliHelper.hasOption("sg");
        boolean skipPhenotypes = cliHelper.hasOption("sp");
        boolean skipAlleles = cliHelper.hasOption("sa");

        // this will download all the data for allele definitions, phenotypes, and drug guidelines
        if (!skipGuidelines || !skipPhenotypes || !skipAlleles) {
          Path zipFile = downloadDir.resolve(sf_zipFileName);
          if (!skipDownload) {
            System.out.println("Downloading data...");
            // download DPWG data
            FileUtils.copyURLToFile(
                new URL("https://s3.pgkb.org/data/" + sf_zipFileName),
                zipFile.toFile());
            ZipUtils.unzip(zipFile, downloadDir);
          }
          if (Files.exists(zipFile)) {
            ZipUtils.unzip(zipFile, downloadDir);
          } else {
            System.out.println("WARNING: Cannot find " + zipFile + " - will have to rely on unpacked content");
          }
        }

        // must get guidelines before alleles
        PgkbGuidelineCollection pgkbGuidelineCollection;
        if (!skipGuidelines) {
          Path drugsDir = cliHelper.getValidDirectory("g", true);
          Path updatedGuidancePath = manager.transformGuidelines(downloadDir, drugsDir);
          pgkbGuidelineCollection = new PgkbGuidelineCollection(updatedGuidancePath);
        } else {
          // if we're skipping new drug data, then use the default data
          pgkbGuidelineCollection = new PgkbGuidelineCollection();
        }

        PhenotypeMap phenotypeMap;
        if (!skipPhenotypes) {
          // transform phenotypes
          Path phenoDir = cliHelper.getValidDirectory("p", true);
          if (!phenoDir.getFileName().endsWith("phenotype")) {
            phenoDir = phenoDir.resolve("phenotype");
          }
          manager.transformPhenotypes(downloadDir, phenoDir);
          phenotypeMap = new PhenotypeMap(phenoDir);
        } else {
          // if we're skipping new phenotype data, then use the default data
          phenotypeMap = new PhenotypeMap();
        }
        validatePhenotypes(phenotypeMap);

        DefinitionReader definitionReader;
        if (!skipAlleles) {
          Path allelesDir = cliHelper.getValidDirectory("a", true);

          Path exemptionsTsv = downloadDir.resolve("exemptions.tsv");
          if (!skipDownload) {
            // download exemptions
            FileUtils.copyURLToFile(new URL(String.format(sf_googleDocUrlFmt, "1xHvvXQIMv3xbqNhuN7zG6WP4DB7lpQDmLvz18w-u_lk")),
                exemptionsTsv.toFile());
          }
          // transform exemptions
          Path exemptionsJson = allelesDir.resolve(EXEMPTIONS_JSON_FILE_NAME);
          Map<String, DefinitionExemption> exemptionsMap = manager.transformExemptions(exemptionsTsv, exemptionsJson);

          // transform allele definitions
          definitionReader = manager.transformAlleleDefinitions(downloadDir, allelesDir, exemptionsMap);

        } else {
          // if we're skipping new gene data, then use the default data
          definitionReader = DefinitionReader.defaultReader();
        }

        List<String> genesUsedInDrugRecommendations = new ArrayList<>(pgkbGuidelineCollection.getGenesWithRecommendations());
        genesUsedInDrugRecommendations.removeAll(definitionReader.getGeneAlleleCount().keySet());
        genesUsedInDrugRecommendations.stream()
            .filter(g -> !g.startsWith("HLA"))
            .map(g -> "WARNING: Gene used in drug recommendation has no allele mapping: " + g)
            .forEach(System.out::println);


        if (cliHelper.hasOption("doc")) {
          Path docsDir = cliHelper.getValidDirectory("doc", true);
          new GeneDrugSummary(definitionReader, phenotypeMap, pgkbGuidelineCollection).write(docsDir);
        }

      } finally {
        if (downloadDir != null && !cliHelper.hasOption("dl")) {
          FileUtils.deleteQuietly(downloadDir.toFile());
        }
      }

    } catch (Exception ex) {
      //noinspection CallToPrintStackTrace
      ex.printStackTrace();
      System.exit(1);
    }
  }

  private Path transformGuidelines(Path downloadDir, Path guidelinesDir) throws IOException {
    if (!Files.exists(guidelinesDir)) {
      Files.createDirectories(guidelinesDir);
    }

    Path downloadedGuidanceFile = downloadDir.resolve(PgkbGuidelineCollection.PRESCRIBING_GUIDANCE_FILE_NAME);
    Path destinationGuidanceFile = guidelinesDir.resolve(PgkbGuidelineCollection.PRESCRIBING_GUIDANCE_FILE_NAME);
    System.out.println();
    System.out.println("Saving guidelines to " + destinationGuidanceFile);

    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection(downloadedGuidanceFile);
    pgkbGuidelineCollection.serializeToJson(destinationGuidanceFile);
    return destinationGuidanceFile;
  }


  private DefinitionFile[] parseDefinitionFiles(Path downloadDir) throws IOException {
    Path definitionsFile = downloadDir.resolve(ALLELES_FILE_NAME);
    if (!Files.exists(definitionsFile)) {
      throw new IOException("Cannot find alleles definitions (" + definitionsFile + ")");
    }
    String json = FileUtils.readFileToString(definitionsFile.toFile(), Charsets.UTF_8);
    return DataSerializer.GSON.fromJson(json, DefinitionFile[].class);
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private DefinitionReader transformAlleleDefinitions(Path downloadDir, Path definitionsDir,
      Map<String, DefinitionExemption> exemptionsMap) throws Exception {

    System.out.println();
    System.out.println("Generating allele definitions...");
    List<DefinitionFile> definitionFiles = new ArrayList<>();
    for (DefinitionFile df : parseDefinitionFiles(downloadDir)) {
      df.setSource(DataSource.PHARMGKB);
      df.validateAlleleNames();
      // strip structural variants since they are unmatchable
      // do this immediately because everything else assumes no structural variants
      InternalWrapper.removeStructuralVariants(df);
      definitionFiles.add(df);
    }

    SortedMap<String, DefinitionFile> definitionFileMap = new TreeMap<>();
    try (VcfHelper vcfHelper = new VcfHelper()) {
      for (DefinitionFile df : definitionFiles) {
        String gene = df.getGeneSymbol();
        DefinitionExemption exemption = exemptionsMap.get(gene);
        if (exemption != null) {
          if (!exemption.getIgnoredAlleles().isEmpty()) {
            System.out.println("Removing ignored named alleles in " + gene + "...");
            InternalWrapper.removeIgnoredNamedAlleles(df, exemption);
          }
          if (!exemption.getIgnoredPositions().isEmpty()) {
            System.out.println("Removing ignored positions in " + gene + "...");
            InternalWrapper.removeIgnoredPositions(df, exemption);
          }
        }

        InternalWrapper.doVcfTranslation(df, vcfHelper);
        definitionFileMap.put(gene, df);
      }
    }

    fixCyp2c19(definitionFileMap.get("CYP2C19"));

    System.out.println();
    System.out.println("Saving allele definitions to " + definitionsDir.toString());
    Set<String> currentFiles = getCurrentFiles(definitionsDir, "_translation.json");

    for (String gene : definitionFileMap.keySet()) {
      DefinitionFile definitionFile = definitionFileMap.get(gene);
      // output file
      Path jsonFile = definitionsDir.resolve(gene + "_translation.json");
      DataSerializer.serializeToJson(definitionFile, jsonFile);
      if (m_verbose) {
        System.out.println("Wrote " + jsonFile);
      }
      if (!currentFiles.remove(gene + "_translation.json")) {
        System.out.println("New gene: " + gene);
      }
    }

    deleteObsoleteFiles(definitionsDir, currentFiles);
    exportVcfData(definitionsDir);

    return new DefinitionReader(definitionsDir);
  }


  private Set<String> getCurrentFiles(Path dir, String suffix) throws IOException {
    Set<String> currentFiles = new HashSet<>();
    try (Stream<Path> list = Files.list(dir)) {
      list.map(PathUtils::getFilename)
          .filter(f -> f.endsWith(suffix))
          .forEachOrdered(currentFiles::add);
    }
    return currentFiles;
  }


  /**
   * Copy any missing alleles from *1 from *38.
   */
  public static void fixCyp2c19(DefinitionFile definitionFile) {
    Preconditions.checkNotNull(definitionFile);
    NamedAllele star1 = Objects.requireNonNull(definitionFile.getNamedAllele("*1"));
    NamedAllele star38 = Objects.requireNonNull(definitionFile.getNamedAllele("*38"));
    star1.initialize(definitionFile.getVariants());
    star38.initialize(definitionFile.getVariants());
    for (int x = 0; x < star38.getAlleles().length; x += 1) {
      if (star1.getAlleles()[x] == null) {
        star1.getAlleles()[x] = star38.getAlleles()[x];
      }
    }
  }


  public static void exportVcfData(Path definitionsDir) throws IOException {
    DefinitionReader definitionReader = new DefinitionReader(definitionsDir);
    SortedSet<String> genes = new TreeSet<>(definitionReader.getGeneAlleleCount().keySet());
    Path positionsFile = definitionsDir.resolve(POSITIONS_VCF);
    System.out.println();
    System.out.println("Saving positions VCF to " + positionsFile);
    VcfHelper.extractPositions(genes, definitionReader, positionsFile);
    Path bgzFile = DockerRunner.bgzip(positionsFile);
    System.out.println("Saved bgzip'd positions VCF to " + bgzFile);
    DockerRunner.prepPharmcatPositions(bgzFile);
  }


  private Map<String, DefinitionExemption> transformExemptions(Path tsvFile, Path jsonFile) throws IOException {
    System.out.println();
    System.out.println("Saving exemptions to " + jsonFile.toString());
    Set<DefinitionExemption> exemptions = m_dataSerializer.deserializeExemptionsFromTsv(tsvFile);
    DataSerializer.serializeToJson(exemptions, jsonFile);

    Map<String, DefinitionExemption> exemptionsMap = new HashMap<>();
    exemptions.forEach(exemption -> exemptionsMap.put(exemption.getGene(), exemption));
    return exemptionsMap;
  }

  private void transformMessages(Path tsvFile, Path jsonFile) throws IOException {
    System.out.println("Saving messages to " + jsonFile.toString());
    DataSerializer.serializeToJson(m_dataSerializer.deserializeMessagesFromTsv(tsvFile), jsonFile);
  }


  private void deleteObsoleteFiles(Path dir, Set<String> obsoleteFilenames) {

    if (!obsoleteFilenames.isEmpty()) {
      System.out.println();
      for (String filename : obsoleteFilenames) {
        Path file = dir.resolve(filename);
        System.out.println("*** Deleting obsolete file: " + file);
        FileUtils.deleteQuietly(file.toFile());
      }
      System.out.println();
    }
  }


  private void transformPhenotypes(Path downloadDir, Path phenoDir) throws IOException {
    Path cpicDir = phenoDir.resolve("cpic");
    System.out.println("Saving CPIC phenotypes to " + cpicDir);
    doTransformPhenotypes(downloadDir.resolve("cpic_phenotypes.json"), cpicDir, DataSource.CPIC);

    Path dpwgDir = phenoDir.resolve("dpwg");
    System.out.println("Saving DPWG phenotypes to " + dpwgDir);
    doTransformPhenotypes(downloadDir.resolve("dpwg_phenotypes.json"), dpwgDir, DataSource.DPWG);
  }

  private void doTransformPhenotypes(Path phenotypeFile, Path outputDir, DataSource source) throws IOException {
    if (!Files.exists(outputDir)) {
      Files.createDirectories(outputDir);
    }

    Set<String> currentFiles = getCurrentFiles(outputDir, ".json");
    try (BufferedReader reader = Files.newBufferedReader(phenotypeFile)) {
      GenePhenotype[] rez = DataSerializer.GSON.fromJson(reader, GenePhenotype[].class);
      Set<String> genes = new HashSet<>();
      for (GenePhenotype gp : rez) {
        if (!genes.add(gp.getGene())) {
          throw new IllegalStateException("Multiple " + source + " GenePhenotypes for " + gp.getGene());
        }

        // generate diplotype data
        gp.generateDiplotypes(source);

        String filename = sanitizeFilename(gp.getGene()) + ".json";
        try (Writer writer = Files.newBufferedWriter(outputDir.resolve(filename))) {
          DataSerializer.GSON.toJson(gp, writer);
        }
        if (!currentFiles.remove(filename)) {
          System.out.println("New gene: " + gp.getGene());
        }
      }
      System.out.println("Found " + rez.length + " " + source + " phenotypes");
      deleteObsoleteFiles(outputDir, currentFiles);
    }
  }


  private static void validatePhenotypes(PhenotypeMap phenotypeMap) {
    for (GenePhenotype gp : phenotypeMap.getCpicGenes()) {
      checkForDuplicatePhenotypeKeys(gp, DataSource.CPIC);
    }
    for (GenePhenotype gp : phenotypeMap.getDpwgGenes()) {
      checkForDuplicatePhenotypeKeys(gp, DataSource.DPWG);
    }
    // validate DPYD phenotypes (DpydCaller depends on this expectation)
    GenePhenotype dpwgGp = Objects.requireNonNull(phenotypeMap.getPhenotype("DPYD", DataSource.DPWG));
    GenePhenotype cpicGp = Objects.requireNonNull(phenotypeMap.getPhenotype("DPYD", DataSource.CPIC));
    for (String hap : dpwgGp.getHaplotypes().keySet()) {
      if (!cpicGp.getHaplotypes().containsKey(hap)) {
        throw new IllegalStateException("DPWG has DPYD " + hap + " but CPIC does not");
      }
    }
  }

  private static void checkForDuplicatePhenotypeKeys(GenePhenotype gp, DataSource source) {
    Set<String> keys = new HashSet<>();
    for (DiplotypeRecord dr : gp.getDiplotypes()) {
      String key = dr.getDiplotypeKey().keySet().stream()
          .sorted()
          .map(h -> h + " (" + dr.getDiplotypeKey().get(h) + ")")
          .collect(Collectors.joining("/"));
      if (!keys.add(key)) {
        throw new IllegalStateException("Duplicate key: " + key + " for " + gp.getGene() + " from " + source);
      }
    }
  }


  private static String sanitizeFilename(String basename) {
    return basename.replaceAll("\\p{Punct}", " ")
        .replaceAll("\\s+", "_");
  }
}
