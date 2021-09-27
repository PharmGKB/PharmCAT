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
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import com.google.common.base.Charsets;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableSet;
import com.google.gson.Gson;
import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.Iupac;
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
  private static final String PHENOTYPES_JSON_FILE_NAME = "gene_phenotypes.json";
  private static final String SUMMARY_REPORT = "summary.md";
  private static final Set<String> PREFER_OUTSIDE_CALL = ImmutableSet.of("CYP2D6", "G6PD", "MT-RNR1");
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
          .addOption("a", "alleles-dir", "directory to save generated allele definition files", true, "a")
          .addOption("m", "messages-dir", "directory to write messages to", true, "m")
          .addOption("g", "drugs-dir", "directory to save drug data to", false, "g")
          .addOption("d", "documenation-dir", "directory to save documentation to", false, "documentation-path")
          .addOption("p", "phenotypes-dir", "directory to save phenotypes to", false, "phenotypes-dir-path")
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

        Map<String,Integer> geneAlleleCountMap;
        if (!skipAlleles) {
          Map<String, DefinitionExemption> exemptionsMap = manager.transformExemptions(exemptionsTsv, exemptionsJson);
          // if we're loading new gene data, use the list of genes from the new data
          geneAlleleCountMap = manager.transformAlleleDefinitions(downloadDir, allelesDir, exemptionsMap);
        } else {
          // if we're skipping new gene data, then use the existing list of genes
          geneAlleleCountMap = new DefinitionReader().getGeneAlleleCount();
        }

        if (!cliHelper.hasOption("sm")) {
          manager.transformMessages(messagesTsv, messagesJson);
        }

        PhenotypeMap phenotypeMap;
        if (cliHelper.hasOption("p")) {
          Path phenoDir = cliHelper.getValidDirectory("p", true);
          phenotypeMap = new PhenotypeMap(manager.writePhenotypes(downloadDir, phenoDir));
        } else {
          phenotypeMap = new PhenotypeMap();
        }

        if (cliHelper.hasOption("d")) {
          Path docsDir = cliHelper.getValidDirectory("d", true);
          manager.writeSummary(docsDir, geneAlleleCountMap, drugs, phenotypeMap);
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

    FileUtils.copyURLToFile(
        new URL("https://files.cpicpgx.org/data/report/current/gene_phenotypes.json"),
        downloadDir.resolve(PHENOTYPES_JSON_FILE_NAME).toFile());

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
    System.out.println("Saving drugs to " + drugsPath);

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
  private Map<String,Integer> transformAlleleDefinitions(Path downloadDir, Path definitionsDir,
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
      if (df.getGeneSymbol().equals("MT-RNR1")) {
        continue;
      }
      doVcfTranslation(df, true);
      definitionFileMap.put(df.getGeneSymbol(), df);
    }

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

    Map<String,Integer> geneAlleleCountMap = new TreeMap<>();
    // output file
    for (String gene : definitionFileMap.keySet()) {
      DefinitionFile definitionFile = definitionFileMap.get(gene);
      geneAlleleCountMap.put(gene, definitionFile.getNamedAlleles().size());
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
    return geneAlleleCountMap;
  }


  private static final Pattern sf_hgvsSnpPattern = Pattern.compile("g\\.\\d+[ACGT]>([ACGT])$");
  private static final Pattern sf_hgvsRepeatPattern = Pattern.compile("g\\.[\\d_]+([ACGT]+\\[\\d+])$");
  private static final Pattern sf_hgvsInsPattern = Pattern.compile("g\\.[\\d_]+(del[ACGT]*)?(ins[ACGT]+)$");
  private static final Pattern sf_hgvsDelPattern = Pattern.compile("g\\.[\\d_]+del[ACGT]*$");

  private void doVcfTranslation(DefinitionFile df, boolean strict) throws IOException {

    NamedAllele referenceNamedAllele = df.getNamedAlleles().stream()
        .filter(NamedAllele::isReference)
        .findFirst()
        .orElseThrow(() -> new IllegalStateException(df.getGeneSymbol() + " does not have reference named allele"));
    referenceNamedAllele.initializeCpicData(df.getVariants());

    VcfHelper vcfHelper = new VcfHelper();
    for (int x = 0; x < df.getVariants().length; x += 1) {
      VariantLocus vl = df.getVariants()[x];
      String errorLocation = df.getGeneSymbol() + " @ " + vl.getCpicPosition();

      String refAllele = Objects.requireNonNull(referenceNamedAllele.getCpicAllele(vl));
      if (isWobble(refAllele)) {
        throw new IllegalStateException(df.getGeneSymbol() + " reference (" + referenceNamedAllele.getName() +
            ") @ " + vl.getCpicPosition() + " is ambiguous: " + refAllele);
      }
      SortedSet<String> altAlleles = new TreeSet<>();
      boolean isSnp = true;
      int numRepeats = 0;
      for (String allele : vl.getCpicAlleles()) {
        if (allele.length() > 1) {
          isSnp = false;
        }
        if (allele.contains("(") || allele.contains(")")) {
          if (!allele.contains("(") || !allele.contains(")")) {
            throw new IllegalStateException(errorLocation + ": allele has mismatched parentheses - " + allele);
          }
          numRepeats += 1;
        }
        if (allele.contains("[") || allele.contains("]")) {
          throw new IllegalStateException(errorLocation + ": allele uses square brackets - " + allele);
        }
        if (!allele.equals(refAllele) && !isWobble(allele)) {
          altAlleles.add(allele);
        }
      }
      if (numRepeats > 0 && numRepeats != vl.getCpicAlleles().size()) {
        throw new IllegalStateException(errorLocation + ": has " + numRepeats + " repeat alleles but " +
            vl.getCpicAlleles().size() + " total alleles (" + vl.getCpicAlleles() + ")");
      }

      List<String> hgvsNames = sf_semicolonSplitter.splitToList(vl.getChromosomeHgvsName());

      if (!isSnp && numRepeats == 0 && altAlleles.size() != 1) {
        // in/dels - must have HGVS to represent each change
        throw new IllegalStateException(errorLocation + ": has " + altAlleles.size() + " alt alleles; max is 1");
      }


      Map<String, String> vcfMap = new HashMap<>();
      List<String> missingAlts = new ArrayList<>(altAlleles);
      long vcfPosition = -1;

      if (isSnp) {
        for (String h : hgvsNames) {
          String hgvs = df.getRefSeqChromosome() + ":" + h;
          VcfHelper.VcfData vcf;
          if (strict) {
            vcf = vcfHelper.hgvsToVcf(hgvs);
          } else {
            Matcher m = sf_hgvsSnpPattern.matcher(h);
            if (!m.matches()) {
              throw new IllegalStateException(errorLocation + ": Invalid substitution HGVS - " + h);
            }
            vcf = new VcfHelper.VcfData(vl.getCpicPosition(), refAllele, m.group(1));
          }

          if (vcfPosition == -1) {
            vcfPosition = vcf.pos;
          } else if (vcfPosition != vcf.pos) {
            throw new IllegalStateException(errorLocation + ": SNP position mismatch (" + vcfPosition + " vs. " +
                vcf.pos + " for " + hgvs + ")");
          }

          if (!refAllele.equals(vcf.ref)) {
            throw new IllegalStateException(errorLocation + ": VCF's reference allele does not match (" + vcfPosition +
                " vs. " + vcf.pos + " for " + hgvs + ")");
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

      } else if (numRepeats > 0) {
        Map<String, VcfHelper.VcfData> firstPass = new HashMap<>();
        for (String h : hgvsNames) {
          Matcher m = sf_hgvsRepeatPattern.matcher(h);
          if (!m.matches()) {
            throw new IllegalStateException(errorLocation + ": Invalid HGVS repeat (" + h + ")");
          }
          String repeatAlt = m.group(1).replaceAll("\\[", "(").replaceAll("]", ")");
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

  private boolean isWobble(String allele) {
    if (allele.length() == 1) {
      return Objects.requireNonNull(Iupac.lookup(allele)).isAmbiguity();
    }
    return false;
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

  private void writeSummary(Path documenationDir, Map<String,Integer> geneAlleleCount, DrugCollection drugs, PhenotypeMap phenotypeMap) throws IOException {
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
        String matcherGeneList = geneAlleleCount.keySet().stream()
            .filter(g -> !PREFER_OUTSIDE_CALL.contains(g))
            .sorted().map(g -> "- " + g + " (" + geneAlleleCount.get(g) + " alleles)")
            .collect(Collectors.joining("\n"));
        String outsideGeneList = PREFER_OUTSIDE_CALL.stream()
            .map(phenotypeMap::lookup)
            .filter(Optional::isPresent)
            .map(Optional::get)
            .map(g -> "- " + g.getGene() + " (" + g.getHaplotypes().size() + " alleles)")
            .collect(Collectors.joining("\n"));
        String drugList = drugs.listReportable().stream().map(d -> "- " + d.getDrugName()).collect(Collectors.joining("\n"));
        IOUtils.write(String.format(mdTemplate, matcherGeneList, outsideGeneList, drugList), fw);
      }
      System.out.println("Saving summary file to " + summaryFile);
    }
  }

  private Path writePhenotypes(Path downloadDir, Path outputDir) throws IOException {
    Path outputPath = outputDir.resolve(PHENOTYPES_JSON_FILE_NAME);
    FileUtils.copyFile(
        downloadDir.resolve(PHENOTYPES_JSON_FILE_NAME).toFile(),
        outputPath.toFile()
    );
    System.out.println("Saving phenotypes to " + outputPath);
    return outputPath;
  }
}
