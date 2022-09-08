package org.pharmgkb.pharmcat.util;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;
import java.util.zip.ZipFile;
import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.gson.Gson;
import org.apache.commons.io.FileUtils;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.MessageList;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.Iupac;
import org.pharmgkb.pharmcat.reporter.DrugCollection;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class manages external resources (e.g. allele definition files, dosing guideline annotations).
 *
 * @author Mark Woon
 */
public class DataManager {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  public static final Path DEFAULT_DEFINITION_DIR = PathUtils.getPathToResource("org/pharmgkb/pharmcat/definition/alleles");
  public static final String EXEMPTIONS_JSON_FILE_NAME = "exemptions.json";
  private static final String PGKB_GUIDELINE_ZIP_FILE_NAME = "guidelineAnnotations.extended.json.zip";
  private static final String POSITIONS_VCF = "pharmcat_positions.vcf";
  private static final String DPWG_ALLELES_FILE_NAME = "dpwg_allele_translations.json";
  private static final String CPIC_ALLELES_FILE_NAME = "allele_definitions.json";
  private static final String sf_dpwgPhenotypesZipFileName = "dpwg_phenotypes.zip";
  private static final String sf_googleDocUrlFmt = "https://docs.google.com/spreadsheets/d/%s/export?format=tsv";

  private static final Splitter sf_semicolonSplitter = Splitter.on(";").trimResults().omitEmptyStrings();
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();
  private final DataSerializer m_dataSerializer = new DataSerializer();
  private final boolean m_verbose;


  private DataManager(boolean verbose) {
    m_verbose = verbose;
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("dl", "download-dir", "directory to save downloaded files", false, "dl")
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

      Path downloadDir;
      if (cliHelper.hasOption("dl")) {
        downloadDir = cliHelper.getValidDirectory("dl", true);
      } else {
        downloadDir = Files.createTempDirectory("pharmcat");
        downloadDir.toFile().deleteOnExit();
      }

      try {
        DataManager manager = new DataManager(cliHelper.isVerbose());

        if (!cliHelper.hasOption("sm")) {
          Path messagesTsv = downloadDir.resolve("messages.tsv");
          // download messages
          FileUtils.copyURLToFile(new URL(String.format(sf_googleDocUrlFmt, "1MkWV6TlJTnw-KRNWeylyUJAUocCgupcJLlmV2fRdtcM")),
              messagesTsv.toFile());

          Path messageDir = cliHelper.getValidDirectory("m", true);
          Path messagesJson = messageDir.resolve(MessageList.MESSAGES_JSON_FILE_NAME);
          manager.transformMessages(messagesTsv, messagesJson);
        }

        // must get guidelines before alleles
        DrugCollection drugs;
        PgkbGuidelineCollection pgkbGuidelineCollection;
        if (!cliHelper.hasOption("sg")) {
          Path drugsDir = cliHelper.getValidDirectory("g", true);
          // download CPIC guidelines
          drugs = DrugCollection.download(drugsDir);
          DrugCollection existingDrugs = new DrugCollection();
          existingDrugs.diff(drugs).forEach(sf_logger::info);

          Path dpwgGuidelinesZip = downloadDir.resolve(PGKB_GUIDELINE_ZIP_FILE_NAME);
          // download DPWG guidelines
          FileUtils.copyURLToFile(
              new URL("https://s3.pgkb.org/data/guidelineAnnotations.extended.json.zip"),
              dpwgGuidelinesZip.toFile());
          // transform DPWG guidelines
          Path dpwgGuidelinesDir = drugsDir.resolve("guidelines");
          manager.transformPgkbGuidelines(dpwgGuidelinesZip, dpwgGuidelinesDir);
          pgkbGuidelineCollection = new PgkbGuidelineCollection(dpwgGuidelinesDir);
        } else {
          // if we're skipping new drug data, then use the default data
          drugs = new DrugCollection();
          pgkbGuidelineCollection = new PgkbGuidelineCollection();
        }

        PhenotypeMap phenotypeMap;
        if (!cliHelper.hasOption("sp")) {
          // download phenotypes
          FileUtils.copyURLToFile(
              new URL("https://files.cpicpgx.org/data/report/current/gene_phenotypes.json"),
              downloadDir.resolve(PhenotypeMap.CPIC_PHENOTYPES_JSON_FILE_NAME).toFile());
          FileUtils.copyURLToFile(
              new URL("https://s3.pgkb.org/data/dpwg_phenotypes.zip"),
              downloadDir.resolve(sf_dpwgPhenotypesZipFileName).toFile());

          // transform phenotypes
          Path phenoDir = cliHelper.getValidDirectory("p", true);
          manager.transformPhenotypes(downloadDir, phenoDir);
          phenotypeMap = new PhenotypeMap(phenoDir);
        } else {
          // if we're skipping new phenotype data, then use the default data
          phenotypeMap = new PhenotypeMap();
        }

        Map<String,Integer> geneAlleleCountMap;
        if (!cliHelper.hasOption("sa")) {
          Path allelesDir = cliHelper.getValidDirectory("a", true);

          Path exemptionsTsv = downloadDir.resolve("exemptions.tsv");
          // download exemptions
          FileUtils.copyURLToFile(new URL(String.format(sf_googleDocUrlFmt, "1xHvvXQIMv3xbqNhuN7zG6WP4DB7lpQDmLvz18w-u_lk")),
              exemptionsTsv.toFile());
          // transform exemptions
          Path exemptionsJson = allelesDir.resolve(EXEMPTIONS_JSON_FILE_NAME);
          Map<String, DefinitionExemption> exemptionsMap = manager.transformExemptions(exemptionsTsv, exemptionsJson);

          // download allele definitions
          FileUtils.copyURLToFile(
              new URL("https://files.cpicpgx.org/data/report/current/allele_definitions.json"),
              downloadDir.resolve(CPIC_ALLELES_FILE_NAME).toFile());
          FileUtils.copyURLToFile(
              new URL("https://s3.pgkb.org/data/dpwg_allele_translations.json"),
              downloadDir.resolve(DPWG_ALLELES_FILE_NAME).toFile());
          // transform allele definitions
          geneAlleleCountMap = manager.transformAlleleDefinitions(downloadDir, allelesDir, exemptionsMap);

        } else {
          // if we're skipping new gene data, then use the default data
          DefinitionReader definitionReader = new DefinitionReader();
          definitionReader.read(DataManager.DEFAULT_DEFINITION_DIR);
          geneAlleleCountMap = definitionReader.getGeneAlleleCount();
        }

        List<String> genesUsedInDrugRecommendations = new ArrayList<>(drugs.list().stream()
            .flatMap(drug -> drug.getGenes().stream())
            .sorted().distinct().toList());
        genesUsedInDrugRecommendations.addAll(pgkbGuidelineCollection.getGenes());
        genesUsedInDrugRecommendations.removeAll(geneAlleleCountMap.keySet());
        genesUsedInDrugRecommendations.stream()
            .filter(g -> !g.startsWith("HLA"))
            .map(g -> "WARNING: Gene used in drug recommendation has no allele mapping: " + g)
            .forEach(System.out::println);


        if (cliHelper.hasOption("doc")) {
          Path docsDir = cliHelper.getValidDirectory("doc", true);
          new GeneDrugSummary().write(docsDir, geneAlleleCountMap, drugs, phenotypeMap);
        }

      } finally {
        if (downloadDir != null) {
          FileUtils.deleteQuietly(downloadDir.toFile());
        }
      }

    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }


  private void transformPgkbGuidelines(Path dpwgGuidelinesZip, Path guidelinesDir) throws IOException {
    if (!Files.exists(guidelinesDir)) {
      Files.createDirectories(guidelinesDir);
    }
    try (ZipFile zipFile = new ZipFile(dpwgGuidelinesZip.toFile())) {
      zipFile.stream().forEach((zipEntry) -> {
        try {
          Path expectedGuideline = PgkbGuidelineCollection.GUIDELINES_DIR.resolve(zipEntry.getName());
          if (Files.exists(expectedGuideline)) {
            Files.copy(
                zipFile.getInputStream(zipEntry),
                guidelinesDir.resolve(zipEntry.getName()),
                StandardCopyOption.REPLACE_EXISTING);
          }
        } catch (IOException ex) {
          throw new RuntimeException("Error copying " + zipEntry, ex);
        }
      });
    }
  }


  private DefinitionFile[] parseDefinitionFiles(Path downloadDir, String fileName) throws IOException {
    Gson gson = new Gson();
    Path definitionsFile = downloadDir.resolve(fileName);
    if (!Files.exists(definitionsFile)) {
      throw new IOException("Cannot find alleles definitions (" + definitionsFile + ")");
    }
    String json = FileUtils.readFileToString(definitionsFile.toFile(), Charsets.UTF_8);
    return gson.fromJson(json, DefinitionFile[].class);
  }


  /**
   * Does the work for stepping through the files and applying the format.
   */
  private Map<String,Integer> transformAlleleDefinitions(Path downloadDir, Path definitionsDir,
      Map<String, DefinitionExemption> exemptionsMap) throws Exception {

    System.out.println("Generating allele definitions...");
    List<DefinitionFile> definitionFiles = new ArrayList<>();
    definitionFiles.addAll(Arrays.asList(parseDefinitionFiles(downloadDir, CPIC_ALLELES_FILE_NAME)));
    definitionFiles.addAll(Arrays.asList(parseDefinitionFiles(downloadDir, DPWG_ALLELES_FILE_NAME)));

    SortedMap<String, DefinitionFile> definitionFileMap = new TreeMap<>();
    try (VcfHelper vcfHelper = new VcfHelper()) {
      for (DefinitionFile df : definitionFiles) {
        // TODO: remove after v1 is released
        if (df.getGeneSymbol().equals("MT-RNR1")) {
          continue;
        }
        doVcfTranslation(df, vcfHelper);
        df.sortPositions();
        definitionFileMap.put(df.getGeneSymbol(), df);
      }
    }

    // transform data as necessary
    System.out.println("Incorporating exemptions...");
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

    System.out.println("Saving allele definitions in " + definitionsDir.toString());
    Set<String> currentFiles = new HashSet<>();
    try (Stream<Path> list = Files.list(definitionsDir)) {
      list.map(PathUtils::getFilename)
          .filter(f -> f.endsWith("_translation.json"))
          .forEachOrdered(currentFiles::add);
    }

    Map<String,Integer> geneAlleleCountMap = new TreeMap<>();
    for (String gene : definitionFileMap.keySet()) {
      DefinitionFile definitionFile = definitionFileMap.get(gene);
      geneAlleleCountMap.put(gene, definitionFile.getNamedAlleles().size());
      // output file
      Path jsonFile = definitionsDir.resolve(gene + "_translation.json");
      m_dataSerializer.serializeToJson(definitionFile, jsonFile);
      if (m_verbose) {
        System.out.println("Wrote " + jsonFile);
      }
      if (!currentFiles.remove(gene + "_translation.json")) {
        System.out.println("New gene: " + gene);
      }
    }

    exportVcfData(definitionsDir);

    deleteObsoleteFiles(definitionsDir, currentFiles);
    return geneAlleleCountMap;
  }


  private static final Pattern sf_hgvsRepeatPattern = Pattern.compile("g\\.[\\d_]+([ACGT]+\\[\\d+])$");
  private static final Pattern sf_hgvsInsPattern = Pattern.compile("g\\.[\\d_]+(del[ACGT]*)?(ins[ACGT]+)$");
  private static final Pattern sf_hgvsDelPattern = Pattern.compile("g\\.[\\d_]+del[ACGT]*$");

  private void doVcfTranslation(DefinitionFile df, VcfHelper vcfHelper) throws IOException {

    NamedAllele referenceNamedAllele = df.getNamedAlleles().stream()
        .filter(NamedAllele::isReference)
        .findFirst()
        .orElseThrow(() -> new IllegalStateException(df.getGeneSymbol() + " does not have reference named allele"));
    referenceNamedAllele.initializeCpicData(df.getVariants());

    for (int x = 0; x < df.getVariants().length; x += 1) {
      VariantLocus vl = df.getVariants()[x];
      String errorLocation = df.getGeneSymbol() + " @ " + vl.getCpicPosition();

      String refAllele = Objects.requireNonNull(referenceNamedAllele.getCpicAllele(vl));
      if (Iupac.isWobble(refAllele)) {
        String hgvs = sf_semicolonSplitter.splitToList(vl.getChromosomeHgvsName()).get(0);
        VcfHelper.VcfData vcf = vcfHelper.hgvsToVcf(df.getRefSeqChromosome() + ":" + hgvs);
        System.out.println(df.getGeneSymbol() + " reference (" + referenceNamedAllele.getName() +
            ") @ " + vl.getCpicPosition() + " is ambiguous (" + refAllele + "):  using " + vcf.ref + " for VCF");
        refAllele = vcf.ref;
      }
      SortedSet<String> altAlleles = new TreeSet<>();
      boolean isSnp = true;
      List<String> repeats = new ArrayList<>();
      List<String> nonRepeats = new ArrayList<>();
      for (String allele : vl.getCpicAlleles()) {
        if (allele.length() > 1) {
          isSnp = false;
        }
        if (allele.contains("(") || allele.contains(")")) {
          if (!allele.contains("(") || !allele.contains(")")) {
            throw new IllegalStateException(errorLocation + ": allele has mismatched parentheses - " + allele);
          }
          repeats.add(allele);
        } else {
          nonRepeats.add(allele);
        }
        if (allele.contains("[") || allele.contains("]")) {
          throw new IllegalStateException(errorLocation + ": allele uses square brackets - " + allele);
        }
        if (!allele.equals(refAllele) && !Iupac.isWobble(allele)) {
          altAlleles.add(allele);
        }
      }
      if (repeats.size() > 0 && repeats.size() != vl.getCpicAlleles().size()) {
        boolean haveSingle = false;
        if (nonRepeats.size() == 1) {
          String repeatedSequence = repeats.get(0);
          repeatedSequence = repeatedSequence.substring(0, repeatedSequence.indexOf("("));
          haveSingle = nonRepeats.get(0).equals(repeatedSequence);
        }
        if (!haveSingle) {
          throw new IllegalStateException(errorLocation + ": has " + repeats.size() + " repeat alleles but " +
              vl.getCpicAlleles().size() + " total alleles (" + vl.getCpicAlleles() + ")");
        }
      }

      List<String> hgvsNames = sf_semicolonSplitter.splitToList(vl.getChromosomeHgvsName());

      if (!isSnp && repeats.size() == 0 && altAlleles.size() != 1) {
        // in/dels - must have HGVS to represent each change
        throw new IllegalStateException(errorLocation + ": has " + altAlleles.size() + " alt alleles; max is 1");
      }


      Map<String, String> vcfMap = new HashMap<>();
      List<String> missingAlts = new ArrayList<>(altAlleles);
      long vcfPosition = -1;

      if (isSnp) {
        for (String h : hgvsNames) {
          String hgvs = df.getRefSeqChromosome() + ":" + h;
          VcfHelper.VcfData vcf = vcfHelper.hgvsToVcf(hgvs);

          if (vcfPosition == -1) {
            vcfPosition = vcf.pos;
          } else if (vcfPosition != vcf.pos) {
            throw new IllegalStateException(errorLocation + ": SNP position mismatch (" + vcfPosition + " vs. " +
                vcf.pos + " for " + hgvs + ")");
          }

          if (!refAllele.equals(vcf.ref)) {
            throw new IllegalStateException(errorLocation + ": VCF's reference allele does not match (" +
                refAllele + " vs. " + vcf.ref + " for " + hgvs + ")");
          }
          if (!missingAlts.remove(vcf.alt)) {
            throw new IllegalStateException(errorLocation + ": VCF's alt allele does not match (expecting " +
                altAlleles + ",  got " + vcf.alt + " for " + hgvs + ")");
          }

          vcfMap.put(refAllele, vcf.ref);
          vcfMap.put(vcf.alt, vcf.alt);
        }
        // warnings
        if (vcfPosition != vl.getCpicPosition()) {
          System.out.println(errorLocation + ": pos/vcf mismatch " + vl.getCpicPosition() +
              " vs. " + vcfPosition);
        }

      } else if (repeats.size() > 0) {
        Map<String, VcfHelper.VcfData> firstPass = new HashMap<>();
        for (String h : hgvsNames) {
          String repeatAlt;
          // treat dups as a form of repeat
          if (h.endsWith("dup")) {
            String repeatedSequence = repeats.get(0);
            repeatedSequence = repeatedSequence.substring(0, repeatedSequence.indexOf("("));
            repeatAlt = repeatedSequence + "(2)";
          } else {
            Matcher m = sf_hgvsRepeatPattern.matcher(h);
            if (!m.matches()) {
              throw new IllegalStateException(errorLocation + ": Invalid HGVS repeat (" + h + ")");
            }
            repeatAlt = m.group(1).replaceAll("\\[", "(").replaceAll("]", ")");
          }
          if (repeatAlt.equals(refAllele)) {
            continue;
          }
          if (!missingAlts.remove(repeatAlt)) {
            throw new IllegalStateException(errorLocation + ": Repeat alt allele does not match (expecting " +
                altAlleles + ",  got " + repeatAlt + ")");
          }

          String hgvs = df.getRefSeqChromosome() + ":" + h;
          VcfHelper.VcfData vcf = vcfHelper.hgvsToVcf(hgvs);
          firstPass.put(repeatAlt, vcf);
        }
        VcfHelper.VcfData vcf = VcfHelper.normalizeRepeats(df.getChromosome(), firstPass.values());

        vcfPosition = vcf.pos;
        vcfMap.put(refAllele, vcf.ref);
        SortedSet<String> repAlts = new TreeSet<>(sf_commaSplitter.splitToList(vcf.alt));
        if (altAlleles.size() != repAlts.size()) {
          throw new IllegalStateException(errorLocation + ": Expected " + altAlleles.size() +
              " repeats but VCF normalization produced " + repAlts.size());
        }
        Iterator<String> repIt = repAlts.iterator();
        for (String alt : altAlleles) {
          vcfMap.put(alt, repIt.next());
        }

      } else {
        // in/del/dup
        for (String h : hgvsNames) {
          String hgvs = df.getRefSeqChromosome() + ":" + h;
          VcfHelper.VcfData vcf = vcfHelper.hgvsToVcf(hgvs);

          String alt;
          Matcher m = sf_hgvsDelPattern.matcher(h);
          if (m.matches()) {
            alt = "del" + refAllele;
          } else if (h.endsWith("dup")) {
            alt = refAllele + refAllele;
          } else if (h.contains("ins")) {
            m = sf_hgvsInsPattern.matcher(h);
            if (!m.matches()) {
              throw new IllegalStateException(errorLocation + ": unsupported ins or delins - " + h);
            }
            alt = m.group(2);
          } else {
            throw new IllegalStateException(errorLocation + ": Unsupported HGVS - " + h);
          }

          vcfPosition = validateVcfPosition(vcfPosition, vcf, errorLocation);
          validateVcfRef(vcfMap, refAllele, vcf, errorLocation);

          if (!missingAlts.remove(alt)) {
            throw new IllegalStateException(errorLocation + ": Alt allele does not match (expecting " + altAlleles +
                ",  got " + alt + ")");
          }
          vcfMap.put(alt, vcf.alt);
        }
      }

      if (missingAlts.size() > 0) {
        throw new IllegalStateException(errorLocation + ": Missing alts " + missingAlts);
      }

      vl.setPosition(vcfPosition);
      vl.setCpicToVcfAlleleMap(vcfMap);
      String vcfRef = vcfMap.get(refAllele);
      vl.setRef(vcfRef);
      altAlleles.forEach((a) -> {
        String vcfAlt = vcfMap.get(a);
        vl.addAlt(vcfAlt);
      });
    }

    for (NamedAllele na : df.getNamedAlleles()) {
      String[] alleles = new String[na.getCpicAlleles().length];
      for (int x = 0; x < alleles.length; x += 1) {
        String cpicAllele = na.getCpicAlleles()[x];
        if (cpicAllele != null) {
          alleles[x] = df.getVariants()[x].getCpicToVcfAlleleMap().get(cpicAllele);
          if (alleles[x] == null) {
            if (Iupac.lookup(cpicAllele).isAmbiguity()) {
              alleles[x] = cpicAllele;
            } else {
              throw new IllegalStateException("Don't know how to translate CPIC allele '" + cpicAllele + "'");
            }
          }
        }
      }
      na.setAlleles(alleles);
    }
  }


  private long validateVcfPosition(long vcfPosition, VcfHelper.VcfData vcf, String errorLocation) {
    if (vcfPosition != -1 && vcfPosition != vcf.pos) {
      throw new IllegalStateException(errorLocation + ": VCF position mismatch (" + vcfPosition + " vs. " + vcf.pos +
          ")");
    }
    return vcf.pos;
  }

  private void validateVcfRef(Map<String, String> vcfMap, String refAllele, VcfHelper.VcfData vcf,
      String errorLocation) {

    if (vcfMap.containsKey(refAllele)) {
      if (!vcf.alt.equals(vcfMap.get(refAllele))) {
        throw new IllegalStateException(errorLocation + ": VCF ref mismatch (" + vcfMap.get(refAllele) + " vs. " +
            vcf.ref);
      }
    } else {
      vcfMap.put(refAllele, vcf.ref);
    }
  }
  

  /**
   * Copy any missing alleles from *1 from *38.
   */
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
    star1.initialize(definitionFile.getVariants());
    star38.initialize(definitionFile.getVariants());
    for (int x = 0; x < star38.getAlleles().length; x += 1) {
      if (star1.getAlleles()[x] == null) {
        star1.getAlleles()[x] = star38.getAlleles()[x];
      }
    }
  }


  private void exportVcfData(Path definitionsDir) throws IOException {

    DefinitionReader definitionReader = new DefinitionReader();
    definitionReader.read(definitionsDir);

    SortedSet<String> genes = new TreeSet<>(definitionReader.getGeneAlleleCount().keySet());

    Path positionsFile = definitionsDir.resolve(POSITIONS_VCF);
    System.out.println("Saving positions VCF to " + positionsFile);
    VcfHelper.extractPositions(genes, definitionReader, positionsFile);
    Path bgzFile = DockerRunner.bgzip(positionsFile);
    System.out.println("Saved bgzip'd positions VCF to " + bgzFile);
    DockerRunner.tabix(bgzFile);
  }


  private Map<String, DefinitionExemption> transformExemptions(Path tsvFile, Path jsonFile) throws IOException {
    System.out.println("Saving exemptions to " + jsonFile.toString());
    Set<DefinitionExemption> exemptions = m_dataSerializer.deserializeExemptionsFromTsv(tsvFile);
    m_dataSerializer.serializeToJson(exemptions, jsonFile);

    Map<String, DefinitionExemption> exemptionsMap = new HashMap<>();
    exemptions.forEach(exemption -> exemptionsMap.put(exemption.getGene(), exemption));
    return exemptionsMap;
  }

  private void transformMessages(Path tsvFile, Path jsonFile) throws IOException {
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

  private void transformPhenotypes(Path downloadDir, Path outputDir) throws IOException {
    Path cpicFile = outputDir.resolve(PhenotypeMap.CPIC_PHENOTYPES_JSON_FILE_NAME);
    System.out.println("Saving CPIC phenotypes to " + cpicFile);
    FileUtils.copyFile(
        downloadDir.resolve(PhenotypeMap.CPIC_PHENOTYPES_JSON_FILE_NAME).toFile(),
        cpicFile.toFile()
    );

    Path dpwgZipFile = downloadDir.resolve(sf_dpwgPhenotypesZipFileName);
    Path dpwgFile = outputDir.resolve(PhenotypeMap.DPWG_PHENOTYPES_JSON_FILE_NAME);
    System.out.println("Saving DPWG phenotypes to " + dpwgFile);
    AtomicBoolean copied = new AtomicBoolean();
    try (ZipFile zipFile = new ZipFile(dpwgZipFile.toFile())) {
      zipFile.stream().forEach((zipEntry -> {
        if (zipEntry.getName().equals(PhenotypeMap.DPWG_PHENOTYPES_JSON_FILE_NAME))
          try {
            Files.copy(
                zipFile.getInputStream(zipEntry),
                dpwgFile,
                StandardCopyOption.REPLACE_EXISTING);
            copied.set(true);
          } catch (IOException ex) {
            throw new RuntimeException("Error copying " + zipEntry, ex);
          }
      }));
    }
    if (!copied.getPlain()) {
      throw new IOException("Could not find " + PhenotypeMap.DPWG_PHENOTYPES_JSON_FILE_NAME + " in " + dpwgZipFile);
    }
  }
}
