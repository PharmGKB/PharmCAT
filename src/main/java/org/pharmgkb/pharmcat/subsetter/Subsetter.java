package org.pharmgkb.pharmcat.subsetter;

import java.io.BufferedReader;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import com.google.common.base.Preconditions;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.TreeMultimap;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellType;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.ComparisonChain;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.InternalWrapper;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.Iupac;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.DataManager;
import org.pharmgkb.pharmcat.util.DataSerializer;
import org.pharmgkb.pharmcat.util.VcfHelper;


/**
 * Subset definitions.
 *
 * @author Mark Woon
 */
public class Subsetter {
  // inputs
  private final Multimap<String, String> m_allowList = HashMultimap.create();
  private final Multimap<String, FunctionData> m_functionOverrides = TreeMultimap.create();
  private final Map<String, Map<String, NamedAllele>> m_extraDefinitions = new HashMap<>();

  private final SortedMap<String, GeneData> m_geneData = new TreeMap<>();

  private final @Nullable DefinitionReader m_definitionReader;
  private @Nullable PhenotypeMap m_phenotypeMap;



  private Subsetter(@Nullable Path defDir, @Nullable Path phenoDir) throws IOException {
    if (defDir == null && phenoDir == null) {
      throw new IllegalArgumentException("Both definition dir and phenotype dir are null");
    }
    if (defDir != null) {
      m_definitionReader = new DefinitionReader(defDir);
    } else {
      m_definitionReader = null;
    }
    if (phenoDir != null) {
      m_phenotypeMap = new PhenotypeMap(phenoDir);
    } else {
      m_phenotypeMap = null;
    }
  }


  private boolean updateDefinitions() {
    if (m_definitionReader == null) {
      System.out.println("Not updating definitions (no definition dir provided).");
      System.out.println();
      return false;
    }

    boolean updated = false;
    if (!m_allowList.isEmpty()) {
      buildSubset();
      updated = true;
    }

    if (!m_extraDefinitions.isEmpty()) {
      if (m_geneData.isEmpty()) {
        allowAll();
      }

      addExtraDefinitions();
      updated = true;
    }

    return updated;
  }


  /**
   * Adds all allele definitions to allow list.
   */
  private void allowAll() {
    Preconditions.checkState(m_definitionReader != null);

    for (String gene : m_definitionReader.getGenes()) {
      DefinitionFile origDefinition = m_definitionReader.getDefinitionFile(gene);

      origDefinition.getNamedAlleles().stream()
          .map(NamedAllele::getName)
          .forEach(na -> m_allowList.put(gene, na));

      GeneData geneData = new GeneData(origDefinition);
      buildGeneData(geneData, origDefinition);
      m_geneData.put(gene, geneData);
    }
  }


  /**
   * Subsets existing allele definition based on allowlist (fills in m_geneData).
   */
  private void buildSubset() {
    Preconditions.checkState(m_definitionReader != null);

    // make sure each gene has a reference haplotype because we require it
    for (String gene : m_allowList.keySet()) {
      m_allowList.put(gene, m_definitionReader.getReferenceHaplotype(gene).getName());
    }

    for (String gene : m_allowList.keySet()) {
      DefinitionFile origDefinition = m_definitionReader.getDefinitionFile(gene);
      GeneData geneData = new GeneData(origDefinition);
      buildGeneData(geneData, origDefinition);
      m_geneData.put(gene, geneData);
    }
  }


  private void buildGeneData(GeneData geneData, DefinitionFile origDefinition) {

    for (String haplotypeName : m_allowList.get(geneData.gene)) {
      NamedAllele hap = origDefinition.getNamedAllele(haplotypeName);
      if (hap == null) {
        continue;
      }
      geneData.haplotypes.add(hap);
      if (!hap.isReference()) {
        if (geneData.gene.equals("CYP2C19") && haplotypeName.equals("*1")) {
          // special handling for CYP2C19*1 - treat it the same as reference
          continue;
        }
        for (int x = 0; x < hap.getCpicAlleles().length; x += 1) {
          if (hap.getCpicAllele(x) != null) {
            geneData.variants.add(origDefinition.getVariants()[x]);
          }
        }
      }
    }

    // look for wobble-only positions
    List<VariantLocus> wobblePositions = new ArrayList<>();
    for (VariantLocus vl : geneData.variants) {
      if (geneData.haplotypes.stream()
          .filter(na -> !na.isReference() && !(geneData.gene.equals("CYP2C19") && na.getName().equals("*1")))
          .allMatch(na -> {
            String allele = na.getCpicAllele(vl);
            return allele == null || Iupac.isWobble(allele);
          })) {
        wobblePositions.add(vl);
      }
    }
    if (!wobblePositions.isEmpty()) {
      System.out.println(geneData.gene + " has " + wobblePositions.size() + " wobble-only position" +
          (wobblePositions.size() > 1 ? "s" : "") + ": " + wobblePositions.stream()
          .map(VariantLocus::toString).collect(Collectors.joining(", ")) + " - REMOVING!");
      wobblePositions.forEach(geneData.variants::remove);
    }

    for (NamedAllele hap : origDefinition.getNamedAlleles()) {
      if (geneData.haplotypes.contains(hap)) {
        continue;
      }
      int numHits = 0;
      List<VariantLocus> missedVariants = new ArrayList<>(wobblePositions);
      for (int x = 0; x < hap.getCpicAlleles().length; x += 1) {
        if (hap.getCpicAllele(x) != null) {
          if (geneData.variants.contains(origDefinition.getVariants()[x])) {
            numHits += 1;
          } else {
            missedVariants.add(origDefinition.getVariants()[x]);
          }
        }
      }
      if (numHits > 0) {
        if (missedVariants.isEmpty()) {
          geneData.subsetHaplotypes.add(hap);
        } else {
          geneData.overlapHaplotypes.add(hap);
          geneData.missedVariants.addAll(missedVariants);
        }
      }
    }
    geneData.calculateUnused();
  }


  private void addExtraDefinitions() {
    for (String gene : m_extraDefinitions.keySet()) {
      m_geneData.values().stream()
          .filter(gd -> gd.gene.equals(gene))
          .findAny()
          .orElseThrow(() -> new IllegalStateException("No gene data for " + gene))
          .extraHaplotypes.addAll(m_extraDefinitions.get(gene).values());
    }
  }


  /**
   * Writes summary of subset.
   */
  private void writeSummary(Path file) throws IOException {
    Preconditions.checkState(m_definitionReader != null);
    Preconditions.checkState(m_phenotypeMap != null);

    if (!Files.exists(file.getParent())) {
      Files.createDirectories(file.getParent());
    }
    new SummaryWriter(m_phenotypeMap, m_geneData).write(file);
  }


  /**
   * Exports updated definition files.
   */
  private void exportDefinitionFiles(Path dir, Path defsDir) throws IOException {
    Preconditions.checkState(m_definitionReader != null);
    if (defsDir == null) {
      System.out.println("No definition output directory specified.  Skipping export.");
      return;
    }

    System.out.println("Saving allele definitions to " + defsDir);
    if (!Files.exists(defsDir)) {
      Files.createDirectories(defsDir);
    } else {
      deleteObsoleteFiles(defsDir, ".json", m_allowList.keySet(), "definition");
      deleteObsoleteFiles(defsDir, ".vcf", m_allowList.keySet(), "definition");
      deleteObsoleteFiles(defsDir, ".bgz", m_allowList.keySet(), "definition");
      deleteObsoleteFiles(defsDir, ".csi", m_allowList.keySet(), "definition");
    }

    Set<DefinitionExemption> exemptions = new HashSet<>();
    try (VcfHelper vcfHelper = new VcfHelper()) {
      for (GeneData gd : m_geneData.values()) {
        // replace all haplotypes in the definition file with allowed haplotypes
        InternalWrapper.resetNamedAlleles(gd.definitionFile, gd.haplotypes);

        // remove ignored positions
        SortedSet<VariantLocus> ignoredPositions = Arrays.stream(gd.definitionFile.getVariants())
            .filter(vl -> !gd.variants.contains(vl))
            .collect(Collectors.toCollection(TreeSet::new));
        InternalWrapper.removeIgnoredPositions(gd.definitionFile, ignoredPositions, false);

        List<String> extras = new ArrayList<>(m_allowList.get(gd.gene));
        for (NamedAllele na : gd.definitionFile.getNamedAlleles()) {
          if (!m_allowList.get(gd.gene).contains(na.getName())) {
            System.out.println("*** SHOULD NOT BE PRESENT: " + na);
          } else {
            extras.remove(na.getName());
          }
        }
        // add new haplotypes
        for (NamedAllele na : gd.extraHaplotypes) {
          if (!m_allowList.get(gd.gene).contains(na.getName())) {
            System.out.println("*** EXTRA ALLELE NOT ON ALLOW LIST: " + na);
          } else {
            extras.remove(na.getName());
            if (gd.definitionFile.getReferenceNamedAllele().getCpicAlleles().length != na.getCpicAlleles().length) {
              throw new IllegalStateException("Extra allele has " + na.getAlleles().length +
                  " alleles, but reference has " + gd.definitionFile.getReferenceNamedAllele().getAlleles().length);
            }
            InternalWrapper.addNamedAllele(gd.definitionFile, na);
          }
        }

        // warn about missing haplotypes
        for (String extra : extras) {
          System.out.println("WARNING: " + gd.gene + " does not have an allele named \"" + extra + "\" - IGNORING!");
        }

        InternalWrapper.doVcfTranslation(gd.definitionFile, vcfHelper);

        if (gd.gene.equals("CYP2C19")) {
          DataManager.fixCyp2c19(gd.definitionFile);
        }

        // export
        Path jsonFile = defsDir.resolve(gd.gene + "_translation.json");
        DataSerializer.serializeToJson(gd.definitionFile, jsonFile);
        //System.out.println("\tWrote " + jsonFile);

        DefinitionExemption exemption = m_definitionReader.getExemption(gd.gene);
        if (exemption != null) {
          if (exemption.getGene().equals("CYP2C9")) {
            if (!m_geneData.containsKey("CYP4F2") || !m_geneData.containsKey("VKORC1")) {
              // remove extra position only used for warfarin
              exemption.getExtraPositions().stream()
                  .filter(p -> "rs12777823".equals(p.getRsid()))
                  .findAny()
                  .ifPresent(p -> exemption.getExtraPositions().remove(p));
            }
          }
          exemptions.add(exemption);
        }
      }
    }

    // write definitions
    Path exemptionsFile = defsDir.resolve(DataManager.EXEMPTIONS_JSON_FILE_NAME);
    DataSerializer.serializeToJson(exemptions, exemptionsFile);
    // generate positions.vcf
    DataManager.exportVcfData(defsDir);
    Files.move(defsDir.resolve(DataManager.POSITIONS_VCF),
        dir.resolve(DataManager.POSITIONS_VCF), StandardCopyOption.REPLACE_EXISTING);
    Files.move(defsDir.resolve(DataManager.POSITIONS_VCF + ".bgz"),
        dir.resolve(DataManager.POSITIONS_VCF + ".bgz"), StandardCopyOption.REPLACE_EXISTING);
    Files.move(defsDir.resolve(DataManager.POSITIONS_VCF + ".bgz.csi"),
        dir.resolve(DataManager.POSITIONS_VCF + ".bgz.csi"), StandardCopyOption.REPLACE_EXISTING);

    Files.move(defsDir.resolve(DataManager.UNIALLELIC_POSITIONS_VCF + ".bgz"),
        dir.resolve(DataManager.UNIALLELIC_POSITIONS_VCF + ".bgz"), StandardCopyOption.REPLACE_EXISTING);
    Files.move(defsDir.resolve(DataManager.UNIALLELIC_POSITIONS_VCF + ".bgz.csi"),
        dir.resolve(DataManager.UNIALLELIC_POSITIONS_VCF + ".bgz.csi"), StandardCopyOption.REPLACE_EXISTING);
  }

  private void exportPhenotypes(Path dir) throws IOException {
    Preconditions.checkState(m_phenotypeMap != null);
    if (dir == null) {
      System.out.println("No definition output directory specified.  Skipping export.");
      return;
    }

    SortedSet<String> genes = new TreeSet<>();
    genes.addAll(m_allowList.keySet());
    genes.addAll(m_functionOverrides.keySet());
    if (genes.isEmpty()) {
      System.out.println("No phenotype changes.  Skipping phenotype export.");
      return;
    }

    System.out.println("Saving phenotypes to " + dir);
    if (!Files.exists(dir)) {
      Files.createDirectories(dir);
    }
    Path cpicDir = dir.resolve("cpic");
    if (!Files.exists(cpicDir)) {
      Files.createDirectories(cpicDir);
    } else {
      deleteObsoleteFiles(cpicDir, ".json", genes, "CPIC phenotypes");
    }

    SortedSet<String> modified = new TreeSet<>(writePhenotypes(genes, m_phenotypeMap.getGenePhenotypes(), cpicDir,
        DataSource.CPIC));

    /*
    Path dpwgDir = dir.resolve("dpwg");
    if (!Files.exists(dpwgDir)) {
      Files.createDirectories(dpwgDir);
    } else {
      deleteObsoleteFiles(dpwgDir, ".json", genes, "DPWG phenotypes");
    }
    */

    Collection<String> diff = CollectionUtils.subtract(m_functionOverrides.keySet(), modified);
    if (!diff.isEmpty()) {
      System.out.println("WARNING: cannot find phenotypes for " + String.join(", ", diff));
    }
    m_phenotypeMap = new PhenotypeMap(dir);
  }

  private Set<String> writePhenotypes(Collection<String> genes, Collection<GenePhenotype> phenotypes, Path dir,
      DataSource src) throws IOException {
    Set<String> changed = new HashSet<>();
    for (GenePhenotype gp : phenotypes) {
      if (genes.contains(gp.getGene())) {
        boolean updatedFunctions = false;
        if (m_functionOverrides.containsKey(gp.getGene())) {
          GeneData geneData = m_geneData.get(gp.getGene());
          // geneData can be null for outside-call only genes (e.g. CYP2D6)
          if (geneData != null && !geneData.extraHaplotypes.isEmpty()) {
            for (NamedAllele na : geneData.extraHaplotypes) {
              FunctionData fd = m_functionOverrides.get(gp.getGene()).stream()
                  .filter(d -> d.allele.equals(na.getName()))
                  .findAny()
                  .orElse(null);
              if (fd == null) {
                // probably an error, but this is not required
                System.out.println("WARNING: extra definition for " + na.getName() + " but no function was specified");
                continue;
              }
              System.out.println("New named allele: " + gp.getGene() + " " + na.getName());
              updatedFunctions = true;
              gp.addHaplotypeRecord(na.getName(), fd.activityScore, fd.function, null);
            }
          }

          for (FunctionData fd : m_functionOverrides.get(gp.getGene())) {
            if (gp.update(fd.allele, fd.activityScore, fd.function, src)) {
              updatedFunctions = true;
            }
          }
          changed.add(gp.getGene());
        }
        // export
        Path jsonFile = dir.resolve(gp.getGene() + ".json");
        if (updatedFunctions) {
          gp.generateDiplotypes();
        }
        DataSerializer.serializeToJson(gp, jsonFile);
        //System.out.println("\tWrote " + jsonFile);
      }
    }
    return changed;
  }

  private void deleteObsoleteFiles(Path dir, String suffix, Collection<String> genes, String desc) throws IOException {
    try (Stream<Path> fileStream = Files.list(dir)) {
      List<Path> files = fileStream
          .filter(f -> {
            String filename = f.getFileName().toString();
            if (!filename.endsWith(suffix)) {
              return false;
            }
            String gene = filename.substring(0, filename.length() - suffix.length());
            return !genes.isEmpty() && !genes.contains(gene);
          })
          .toList();
      if (!files.isEmpty()) {
        System.out.println("Deleting unused " + desc + ":");
        for (Path file : files) {
          System.out.println("  * " + file.getFileName());
          Files.delete(file);
        }
      }
    }
  }


  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addVersion("PharmCAT Definition Subsetter" + CliUtils.getVersion())
          .addOption("a", "alleles", "Allele whitelist file (tsv or xlsx)", false, "file")
          .addOption("pos", "positions", "Positions whitelist file (bim file)", false, "file")
          .addOption("d", "extra-definitions", "Extra named allele definitions", false, "file")
          .addOption("pc", "cpic-phenotypes", "CPIC phenotype override file (xlsx)", false, "file")
          // input for original/base data - required for subsetting positions
          .addOption("i", "input-data-dir", "Input data directory (base data)", true, "directory")
          .addOption("o", "output-dir", "Output directory (PharmCAT root dir)", false, "directory")
      ;
      if (!cliHelper.parse(args)) {
        if (!cliHelper.hasOption("i")) {
          System.out.println("""
              
              Please specify base data directory with -i.
              The base data directory should have the data from which to generate the subset:
              * copy of the contents of org/pharmgkb/pharmcat/definition/alleles
                to dataDir/definitions.
              * copy org/pharmgkb/pharmcat/phenotype/cpic and org/pharmgkb/pharmcat/phenotype/cpic
                to dataDir/phenotypes.
              """);
        }
        return;
      }

      if (cliHelper.hasOption("pos") && cliHelper.hasOption("a")) {
        System.out.println("-pos and -a are mutually exclusive.");
        return;
      }

      Path dataDir = cliHelper.getValidDirectory("i", true);
      Path baseDefDir = dataDir.resolve("definition");
      if (cliHelper.hasOption("d") || cliHelper.hasOption("pos") || cliHelper.hasOption("a")) {
        if (!Files.isDirectory(baseDefDir)) {
          System.out.println("Cannot find 'definitions' subdirectory in " + dataDir);
        }
      }
      Path basePhenoDir = dataDir.resolve("phenotype");
      if (cliHelper.hasOption("pc")) {
        if (!Files.isDirectory(basePhenoDir)) {
          System.out.println("Cannot find 'phenotypes' subdirectory in " + dataDir);
        }
      }

      Subsetter subsetter = new Subsetter(baseDefDir, basePhenoDir);

      // read in extra named alleles
      if (cliHelper.hasOption("d")) {
        Path file = cliHelper.getValidFile("d", true);
        subsetter.parseExtraDefinitions(file);
      }

      // read in allowList and reportAsReference
      if (cliHelper.hasOption("pos")) {
        Path file = cliHelper.getValidFile("pos", true);
        if (!file.toString().endsWith(".bim")) {
          throw new ReportableException("Unsupported file type for positions whitelist: " + file);
        }
        subsetter.parsePositionsBim(file);

      } else if (cliHelper.hasOption("a")) {
        Path file = cliHelper.getValidFile("a", true);
        if (file.toString().endsWith(".tsv")) {
          subsetter.parseAlleleAllowListTsv(file);
        } else if (file.toString().endsWith(".xlsx")) {
          subsetter.parseAlleleAllowListXls(file);
        } else {
          throw new ReportableException("Unsupported file type for alleles whitelist: " + file);
        }
      }

      // read in phenotype overrides
      if (cliHelper.hasOption("pc")) {
        Path file = cliHelper.getValidFile("pc", true);
        subsetter.parseFunctionOverrides(file);
      }

      if (!subsetter.updateDefinitions()) {
        System.out.println("No definitions updated.");
      }

      if (cliHelper.hasOption("o")) {
        Path outDir = cliHelper.getValidDirectory("o", true);

        Path defDir = outDir.resolve("src/main/resources/org/pharmgkb/pharmcat/definition/alleles");
        if (!Files.isDirectory(defDir)) {
          Files.createDirectories(defDir);
        }
        subsetter.exportDefinitionFiles(outDir, defDir);

        Path phenoDir = outDir.resolve("src/main/resources/org/pharmgkb/pharmcat/phenotype");
        if (!Files.isDirectory(phenoDir)) {
          Files.createDirectories(phenoDir);
        }
        subsetter.exportPhenotypes(phenoDir);

        // write a summary after exporting everything
        subsetter.writeSummary(defDir.resolve("summary.xlsx"));
      }

      System.out.println("Done.");

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      System.out.println(ex.getMessage());
    } catch (Exception e) {
      //noinspection CallToPrintStackTrace
      e.printStackTrace();
    }
  }


  private void parseExtraDefinitions(Path jsonFile) throws IOException {
    Preconditions.checkState(m_definitionReader != null);

    try (BufferedReader reader = Files.newBufferedReader(jsonFile, StandardCharsets.UTF_8)) {
      @SuppressWarnings({ "rawtypes", "unchecked" })
      Map<String, List<Map>> data = DataSerializer.GSON.fromJson(reader, Map.class);
      for (String gene : data.keySet()) {
        @SuppressWarnings({ "rawtypes" })
        List<Map> extraAlleles = data.get(gene);
        for (@SuppressWarnings("rawtypes") Map ea : extraAlleles) {

          @SuppressWarnings("unchecked")
          NamedAllele na = new NamedAllele((String)ea.get("id"), (String)ea.get("name"), null,
              ((List<String>)ea.get("cpicAlleles")).toArray(new String[0]), (boolean)ea.get("reference"));
          m_extraDefinitions.computeIfAbsent(gene, g -> new HashMap<>())
              .put(na.getName(), na);
        }
      }
    }
  }

  private void parsePositionsBim(Path file) throws IOException {
    Preconditions.checkState(m_definitionReader != null);

    Set<String> rsidsOfInterest = m_definitionReader.getLocationsOfInterest().values().stream()
        .map(VariantLocus::getRsid)
        .filter(Objects::nonNull)
        .collect(Collectors.toSet());

    // get allowlisted positions
    Multimap<String, String> foundRsids = HashMultimap.create();
    try (BufferedReader reader = Files.newBufferedReader(file)) {
      String line = reader.readLine();
      while (line != null) {
        String[] data = line.split("\t");
        String rsid = StringUtils.stripToNull(data[1]);
        if (rsid != null && rsidsOfInterest.contains(rsid)) {
          String allele = StringUtils.stripToNull(data[4]);
          if (allele != null) {
            foundRsids.put(rsid, allele);
          }
          allele = StringUtils.stripToNull(data[5]);
          if (allele != null) {
            foundRsids.put(rsid, allele);
          }
        }
        line = reader.readLine();
      }
    }

    // derive allowlisted alleles from positions
    for (String gene : m_definitionReader.getGenes()) {
      VariantLocus[] positions = m_definitionReader.getPositions(gene);

      for (NamedAllele na : m_definitionReader.getHaplotypes(gene)) {
        boolean allMatch = true;
        for (int x = 0; x < positions.length; x += 1) {
          if (na.getCpicAllele(x) != null) {
            VariantLocus vl = positions[x];
            if (!foundRsids.keySet().contains(vl.getRsid()) ||
                !foundRsids.get(vl.getRsid()).contains(na.getCpicAllele(x))) {
              allMatch = false;
              break;
            }
          }
        }
        if (allMatch) {
          m_allowList.put(gene, na.getName());
        }
      }
    }
  }


  private void parseAlleleAllowListTsv(Path tsvFile) throws IOException {
    Preconditions.checkState(m_definitionReader != null);

    try (BufferedReader reader = Files.newBufferedReader(tsvFile)) {
      String line = reader.readLine();
      while (line != null) {
        String[] data = line.split("\t");
        String gene = StringUtils.stripToNull(data[0]);
        if (gene == null || gene.equalsIgnoreCase("gene")) {
          // ignore headers
          line = reader.readLine();
          continue;
        }
        String allele = StringUtils.stripToNull(data[1]);
        if (allele == null) {
          line = reader.readLine();
          continue;
        }
        if (allele.startsWith("\"") && allele.endsWith("\"")) {
          allele = allele.substring(1, allele.length() - 1);
        }
        m_allowList.put(gene, allele);

        if (data.length > 2) {
          String mod = StringUtils.stripToNull(data[2]);
          if (mod != null) {
            System.out.println("Don't know what to do with '" + mod + "' for " + gene + " " + allele);
          }
        }

        line = reader.readLine();
      }
    }
  }

  private void parseAlleleAllowListXls(Path xlsxFile) throws IOException {
    Preconditions.checkState(m_definitionReader != null);

    try (Workbook workbook = WorkbookFactory.create(xlsxFile.toFile())) {
      Sheet sheet = workbook.getSheetAt(0);
      for (Row row : sheet) {
        Cell geneCell = row.getCell(0);
        String gene = StringUtils.stripToNull(geneCell.getStringCellValue());
        if (gene == null || gene.equalsIgnoreCase("gene")) {
          // ignore headers
          continue;
        }

        Cell alleleCell = row.getCell(1);
        String allele = StringUtils.stripToNull(alleleCell.getStringCellValue());
        if (allele == null) {
          continue;
        }
        if (allele.startsWith("\"") && allele.endsWith("\"")) {
          allele = allele.substring(1, allele.length() - 1);
        }

        DefinitionFile definitionFile = m_definitionReader.getDefinitionFile(gene);
        if (definitionFile.getNamedAllele(allele) == null && m_extraDefinitions.get(gene).get(allele) == null) {
          StringBuilder builder = new StringBuilder()
              .append("WARNING: ")
              .append(gene)
              .append(" does not have an allele named \"")
              .append(allele)
              .append("\" listed in the definitions file");
          if (!m_extraDefinitions.isEmpty()) {
            builder.append(" or the extra definitions file");
          }
          System.out.println(builder);
        }
        m_allowList.put(gene, allele);

        Cell modCell = row.getCell(2);
        if (modCell != null) {
        String mod = StringUtils.stripToNull(modCell.getStringCellValue());
          if (mod != null) {
            System.out.println("Don't know what to do with '" + mod + "' for " + gene + " " + allele);
          }
        }
      }
    }
  }


  private void parseFunctionOverrides(Path xlsxFile) throws IOException {
    Preconditions.checkState(m_phenotypeMap != null);

    try (Workbook workbook = WorkbookFactory.create(xlsxFile.toFile())) {
      Sheet sheet = workbook.getSheetAt(0);
      for (Row row : sheet) {
        //System.out.println("Row " + row.getRowNum());
        Cell geneCell = row.getCell(0);
        String gene = StringUtils.stripToNull(geneCell.getStringCellValue());
        if (gene == null || gene.equalsIgnoreCase("gene")) {
          // ignore headers
          continue;
        }

        Cell alleleCell = row.getCell(1);
        String allele = StringUtils.stripToNull(alleleCell.getStringCellValue());
        if (allele == null) {
          continue;
        }
        if (allele.startsWith("\"") && allele.endsWith("\"")) {
          allele = allele.substring(1, allele.length() - 1);
        }

        Cell activityCell = row.getCell(2);
        String activity = null;
        if (activityCell != null) {
          if (activityCell.getCellType() == CellType.NUMERIC) {
            double val = activityCell.getNumericCellValue();
            if (val == 1) {
              activity = "1.0";
            } else if (val == 0) {
              activity = "0.0";
            } else if (val == 2) {
              activity = "2.0";
            } else {
              activity = Double.toString(val);
            }
          } else {
            activity = StringUtils.stripToNull(activityCell.getStringCellValue());
          }
          if (activity != null && activity.startsWith("\"") && activity.endsWith("\"")) {
            activity = activity.substring(1, activity.length() - 1);
          }
        }

        Cell functionCell = row.getCell(3);
        String function = StringUtils.stripToNull(functionCell.getStringCellValue());
        if (function != null && function.startsWith("\"") && function.endsWith("\"")) {
          function = function.substring(1, function.length() - 1);
        }
        m_functionOverrides.put(gene, new FunctionData(gene, allele, activity, function));
      }
    }
  }



  private record FunctionData(String gene, String allele, @Nullable String activityScore, String function)
      implements Comparable<FunctionData> {

    private FunctionData(String gene, String allele, String activityScore, String function) {
      this.gene = gene;
      this.allele = allele;
      this.activityScore = activityScore;
      this.function = function;
    }

    @Override
    public int compareTo(FunctionData o) {
      return new ComparisonChain()
          .compare(gene, o.gene)
          .compare(allele, o.allele)
          .result();
    }
  }
}
