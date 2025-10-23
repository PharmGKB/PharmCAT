package org.pharmgkb.pharmcat.stats;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.lang.invoke.MethodHandles;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.google.common.base.Splitter;
import org.apache.commons.lang3.StringUtils;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.common.util.StreamUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.CombinationMatcher;
import org.pharmgkb.pharmcat.reporter.format.CallsOnlyFormat;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.PoiUtils;

import static org.pharmgkb.pharmcat.reporter.format.CallsOnlyFormat.HEADER_SAMPLE_ID;
import static org.pharmgkb.pharmcat.util.PoiUtils.writeCell;


/**
 * This tool generates allele frequency stats from *.report.tsv files.
 *
 * @author Mark Woon
 */
public class CalcAlleleFrequencies {
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();
  private final Env m_env;
  private final int m_pivotCol;
  private final boolean m_treatCombinationsAsUnknown;
  private final @Nullable String m_gene;
  private final Map<String, GeneStats> m_stats = new TreeMap<>();
  private final NumberFormat m_numberFormat = NumberFormat.getNumberInstance();
  private final NumberFormat m_percentFormat = NumberFormat.getPercentInstance();


  public static void main(@Nullable String[] args) throws Exception {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("i", "input", "Calls-only TSV file or directory", true, "file")
          .addOption("o", "output-dir", "Output to (optional, default is input file directory)", false, "directory")
          .addOption("pc", "pivot-column", "Pivot column number; first column is 1", false, "number")
          .addOption("uc", "unknown-combos", "Treat combination calls as unknown")
          .addOption("g", "gene", "Only generate stats for specified gene", false, "gene")
          ;
      if (!cliHelper.parse(args)) {
        if (!cliHelper.isHelpRequested() && !cliHelper.isVersionRequested()) {
          CliUtils.failIfNotTest();
        }
        return;
      }

      Path input = cliHelper.getPath("i");
      if (!Files.exists(input)) {
        throw new ReportableException("Input file/directory does not exist: " + input);
      }
      Path outDir;
      if (cliHelper.hasOption("o")) {
        outDir = cliHelper.getValidDirectory("o", true);
      } else if (Files.isDirectory(input)) {
        outDir = input;
      } else {
        outDir = input.getParent();
        if (outDir == null) {
          outDir = input.getFileSystem().getPath(".");
        }
      }
      int pivotCol = -1;
      if (cliHelper.hasOption("pc")) {
        // ask for 1-based, but we use 0-based
        pivotCol = cliHelper.getIntValue("pc") - 1;
        if (pivotCol < 16) {
          throw new ReportableException("Invalid data column '" + cliHelper.getIntValue("pc") +
              "' (expecting value > 16)");
        }
      }

      CalcAlleleFrequencies calcAlleleFrequencies = new CalcAlleleFrequencies(new Env(), pivotCol,
          cliHelper.hasOption("uc"), cliHelper.getValue("g"));
      if (Files.isDirectory(input)) {
        calcAlleleFrequencies.ingestDir(input);
      } else if (input.toString().endsWith(".tsv") || input.toString().endsWith(".tsv.gz") ||
          input.toString().endsWith(".tsv.zip")) {
        calcAlleleFrequencies.ingestFile(input);
      } else {
        throw new ReportableException("Input file must be a *.tsv, *.tsv.gz or *.tsv.zip file");
      }
      calcAlleleFrequencies.write(outDir);

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      CliUtils.failIfNotTest(ex.getMessage());
    } catch (Exception ex) {
      CliUtils.failIfNotTest(ex);
    }
  }


  private CalcAlleleFrequencies(Env env, int pivotCol, boolean treatCombinationsAsUnknown, @Nullable String gene) {
    m_env = env;
    m_pivotCol = pivotCol;
    m_treatCombinationsAsUnknown = treatCombinationsAsUnknown;
    m_percentFormat.setMinimumFractionDigits(3);
    m_gene = gene;
  }

  private boolean doPivot() {
    return m_pivotCol >= 0;
  }


  private void ingestDir(Path dir) throws IOException {
    try (DirectoryStream<Path> stream = Files.newDirectoryStream(dir)) {
      for (Path path : stream) {
        if (Files.isDirectory(path)) {
          ingestDir(path);
        } else if (Files.isRegularFile(path) && (path.toString().endsWith(".report.tsv") ||
            path.toString().endsWith(".report.tsv.gz") ||
            path.toString().endsWith(".report.zip"))) {
          ingestFile(path);
        }
      }
    }
  }

  private void ingestFile(Path file) throws IOException {

    System.out.println("Reading " + file);
    int geneCol = 0;
    try (BufferedReader reader = StreamUtils.openReader(file)) {
      String line = reader.readLine();
      int lineNum = 1;
      // skip PharmCAT version info
      if (line.startsWith("PharmCAT")) {
        line = reader.readLine();
        lineNum += 1;
      }
      // check header line
      if (line.startsWith(HEADER_SAMPLE_ID)) {
        geneCol = 1;
      } else if (!line.startsWith("Gene\t")) {
        throw new IOException("Expecting headers but line " + lineNum + " is: " + line);
      }
      List<String> headers = List.of(line.split("\t"));
      int missingPositionsCol = 12;
      if (headers.indexOf(CallsOnlyFormat.HEADER_VARIANTS) > 0 &&
          headers.indexOf(CallsOnlyFormat.HEADER_VARIANTS) < 17) {
        missingPositionsCol = 13;
      }
      if (m_pivotCol >= 0) {
        System.out.println("Pivot column is: " + headers.get(m_pivotCol));
      }

      while ((line = reader.readLine()) != null) {
        lineNum += 1;
        if (StringUtils.isBlank(line)) {
          continue;
        }
        try {
          String[] fields = line.split("\t");
          String gene = fields[geneCol];
          if (m_gene != null && !m_gene.equals(gene)) {
            continue;
          }
          String dip = fields[geneCol + 1];
          String pheno = fields[geneCol + 2];
          String missingPositions = fields[geneCol + missingPositionsCol];
          String bgData = null;
          if (doPivot()) {
            if (m_pivotCol < fields.length) {
              bgData = fields[m_pivotCol];
            }
          }
          GeneStats geneStats = m_stats.computeIfAbsent(gene, GeneStats::new);
          geneStats.add(dip, pheno, bgData, missingPositions);
        } catch (Exception ex) {
          System.out.println("Error on line " + lineNum);
          System.out.println(line);
          throw ex;
        }
      }
    }
  }

  private void write(Path outDir) throws IOException {
    for (String gene : m_stats.keySet()) {
      GeneStats geneStats = m_stats.get(gene);
      Path file = outDir.resolve(gene + "_stats.xlsx");
      try (XSSFWorkbook workbook = new XSSFWorkbook()) {
        CellStyle headerStyle = PoiUtils.getHeaderCellStyle(workbook);
        writeDiplotypeTab(workbook, headerStyle, geneStats);
        writeAlleleTab(workbook, headerStyle, geneStats);
        writePhenotypeTab(workbook, headerStyle, geneStats);
        writeMiscTab(workbook, headerStyle, geneStats);

        System.out.println("Writing to " + file);
        try (OutputStream fos = Files.newOutputStream(file)) {
          workbook.write(fos);
        }
      }
    }
  }


  private void writeDiplotypeTab(Workbook workbook, CellStyle headerStyle, GeneStats geneStats) {
    Sheet sheet = workbook.createSheet("Diplotypes");

    int rowNum = 0;
    Row headerRow = sheet.createRow(rowNum);
    writeCell(headerRow, 0, "Diplotype", headerStyle);
    writeCell(headerRow, 1, "Count (n=" + m_numberFormat.format(geneStats.numSamples) + ")", headerStyle);
    writeCell(headerRow, 2, "Frequency", headerStyle);
    writeBgHeaders(headerRow, headerStyle, geneStats.numSamples, geneStats.numBgSamples, geneStats.diplotypeRegions);

    List<Map.Entry<String, Integer>> entries = new ArrayList<>(geneStats.diplotypeCounts.entrySet());
    // sort entries by values (ascending order)
    entries.sort(Map.Entry.<String, Integer>comparingByValue().reversed());

    for (Map.Entry<String, Integer> entry : entries) {
      String key = entry.getKey();
      rowNum += 1;
      Row row = sheet.createRow(rowNum);

      writeRow(row, geneStats.diplotypeCounts.get(key), geneStats.numSamples, geneStats.numBgSamples,
          geneStats.bgDiplotypeCounts, geneStats.diplotypeRegions, key);
    }
    for (int x = 0; x < headerRow.getLastCellNum(); x += 1) {
      sheet.autoSizeColumn(x);
    }
  }

  private void writeAlleleTab(Workbook workbook, CellStyle headerStyle, GeneStats geneStats) {
    Sheet sheet = workbook.createSheet("Alleles");

    int rowNum = 0;
    Row headerRow = sheet.createRow(rowNum);
    writeCell(headerRow, 0, "Allele", headerStyle);
    writeCell(headerRow, 1, "Count (n=" + m_numberFormat.format(geneStats.numAlleles) + ")", headerStyle);
    writeCell(headerRow, 2, "Frequency", headerStyle);
    writeBgHeaders(headerRow, headerStyle, geneStats.numAlleles, geneStats.numBgAlleles, geneStats.alleleRegions);

    List<Map.Entry<String, Integer>> entries = new ArrayList<>(geneStats.alleleCounts.entrySet());
    // sort entries by values (ascending order)
    entries.sort(Map.Entry.<String, Integer>comparingByValue().reversed());

    for (Map.Entry<String, Integer> entry : entries) {
      String key = entry.getKey();
      rowNum += 1;
      Row row = sheet.createRow(rowNum);
      writeRow(row, geneStats.alleleCounts.get(key), geneStats.numAlleles, geneStats.numBgAlleles,
          geneStats.bgAlleleCounts, geneStats.alleleRegions, key);
    }
    for (int x = 0; x <= 2; x += 1) {
      sheet.autoSizeColumn(x);
    }
  }

  private void writePhenotypeTab(Workbook workbook, CellStyle headerStyle, GeneStats geneStats) {
    Sheet sheet = workbook.createSheet("Phenotypes");

    int rowNum = 0;
    Row headerRow = sheet.createRow(rowNum);
    writeCell(headerRow, 0, "Phenotype", headerStyle);
    writeCell(headerRow, 1, "Count (n=" + m_numberFormat.format(geneStats.numSamples) + ")", headerStyle);
    writeCell(headerRow, 2, "Frequency", headerStyle);
    writeBgHeaders(headerRow, headerStyle, geneStats.numSamples, geneStats.numBgSamples, geneStats.phenotypeRegions);

    List<Map.Entry<String, Integer>> entries = new ArrayList<>(geneStats.phenotypeCounts.entrySet());
    // sort entries by values (ascending order)
    entries.sort(Map.Entry.<String, Integer>comparingByValue().reversed());

    for (Map.Entry<String, Integer> entry : entries) {
      String key = entry.getKey();
      rowNum += 1;
      Row row = sheet.createRow(rowNum);

      writeRow(row, geneStats.phenotypeCounts.get(key), geneStats.numSamples, geneStats.numBgSamples,
          geneStats.bgPhenotypeCounts, geneStats.phenotypeRegions, key);
    }
    for (int x = 0; x < headerRow.getLastCellNum(); x += 1) {
      sheet.autoSizeColumn(x);
    }
  }

  private void writeMiscTab(Workbook workbook, CellStyle headerStyle, GeneStats geneStats) {
    Sheet sheet = workbook.createSheet("Misc");

    int rowNum = 0;
    Row row = sheet.createRow(rowNum);
    writeCell(row, 1, "Total", headerStyle);
    writeCell(row, 2, "Unknown", headerStyle);
    writeCell(row, 3, "Multiple Calls", headerStyle);
    rowNum += 1;

    row = sheet.createRow(rowNum);
    writeCell(row, 0, "Number of samples", headerStyle);
    writeCell(row, 1, m_numberFormat.format(geneStats.numSamples));
    writeCell(row, 2, m_numberFormat.format(geneStats.numUnknownSamples));
    writeCell(row, 3, m_numberFormat.format(geneStats.numMultiDipSamples));
    rowNum += 1;

    row = sheet.createRow(rowNum);
    writeCell(row, 0, "Number of sample with all positions", headerStyle);
    writeCell(row, 1, m_numberFormat.format(geneStats.numSamplesComplete));
    writeCell(row, 2, m_numberFormat.format(geneStats.numUnknownSamplesComplete));
    writeCell(row, 3, m_numberFormat.format(geneStats.numMultiDipSamplesComplete));
    rowNum += 1;

    row = sheet.createRow(rowNum);
    writeCell(row, 0, "Number of samples with missing positions", headerStyle);
    writeCell(row, 1, m_numberFormat.format(geneStats.numSamplesMissingPositions));
    writeCell(row, 2, m_numberFormat.format(geneStats.numUnknownSamplesMissingPositions));
    writeCell(row, 3, m_numberFormat.format(geneStats.numMultiDipSamplesMissingPositions));
    rowNum += 1;

    row = sheet.createRow(rowNum);
    writeCell(row, 0, "Number of samples with undocumented positions", headerStyle);
    writeCell(row, 1, m_numberFormat.format(geneStats.numSamplesUndocumentedVariants));
    writeCell(row, 2, m_numberFormat.format(geneStats.numUnknownSamplesUndocumentedVariants));
    writeCell(row, 3, m_numberFormat.format(geneStats.numMultiDipSamplesUndocumentedVariants));
    rowNum += 1;

    row = sheet.createRow(rowNum);
    writeCell(row, 0, "Number of samples with missing and undocumented positions", headerStyle);
    writeCell(row, 1, m_numberFormat.format(geneStats.numSamplesMissingAndUndocumented));
    writeCell(row, 2, m_numberFormat.format(geneStats.numUnknownSamplesMissingAndUndocumented));
    writeCell(row, 3, m_numberFormat.format(geneStats.numMultiDipSamplesMissingAndUndocumented));
    rowNum += 1;

    if (!geneStats.missingPositions.isEmpty()) {
      rowNum += 1;
      row = sheet.createRow(rowNum);
      writeCell(row, 0, "Missing positions", headerStyle);
      rowNum += 1;
      SortedSet<String> variants = new TreeSet<>(geneStats.missingPositions.keySet());
      variants.addAll(geneStats.unknownMissingPositions.keySet());
      variants.addAll(geneStats.multiDipMissingPositions.keySet());
      for (String variant : variants) {
        row = sheet.createRow(rowNum);
        writeCell(row, 0, variant);
        writeCell(row, 1, m_numberFormat.format(geneStats.missingPositions.getOrDefault(variant, 0)));
        writeCell(row, 2, m_numberFormat.format(geneStats.unknownMissingPositions.getOrDefault(variant, 0)));
        writeCell(row, 3, m_numberFormat.format(geneStats.multiDipMissingPositions.getOrDefault(variant, 0)));
        rowNum += 1;
      }
    }

    if (!geneStats.undocumentedVariants.isEmpty()) {
      rowNum += 1;
      row = sheet.createRow(rowNum);
      writeCell(row, 0, "Undocumented variants", headerStyle);
      rowNum += 1;
      SortedSet<String> variants = new TreeSet<>(geneStats.undocumentedVariants.keySet());
      variants.addAll(geneStats.unknownUndocumentedVariants.keySet());
      variants.addAll(geneStats.multiDipUndocumentedVariants.keySet());
      for (String variant : variants) {
        row = sheet.createRow(rowNum);
        String rsid = geneStats.undocumentedRsids.get(variant);
        if (rsid == null) {
          rsid = "";
        } else {
          rsid += ": ";
        }
        writeCell(row, 0, rsid + variant);
        writeCell(row, 1, m_numberFormat.format(geneStats.undocumentedVariants.getOrDefault(variant, 0)));
        writeCell(row, 2, m_numberFormat.format(geneStats.unknownUndocumentedVariants.getOrDefault(variant, 0)));
        writeCell(row, 3, m_numberFormat.format(geneStats.multiDipUndocumentedVariants.getOrDefault(variant, 0)));
        rowNum += 1;
      }
    }

    if (!geneStats.unknownDipCombinations.isEmpty()) {
      rowNum += 1;
      row = sheet.createRow(rowNum);
      writeCell(row, 0, "Diplotypes converted to Unknown/Unknown", headerStyle);
      rowNum += 1;
      for (String dip : geneStats.unknownDipCombinations.keySet()) {
        row = sheet.createRow(rowNum);
        writeCell(row, 0, dip);
        writeCell(row, 2, m_numberFormat.format(geneStats.unknownDipCombinations.get(dip)));
        rowNum += 1;
      }
    }

    if (!geneStats.multiDipCombinations.isEmpty()) {
      rowNum += 1;
      row = sheet.createRow(rowNum);
      writeCell(row, 0, "Multiple calls", headerStyle);
      rowNum += 1;
      for (String dip : geneStats.multiDipCombinations.keySet()) {
        row = sheet.createRow(rowNum);
        writeCell(row, 0, dip);
        writeCell(row, 3, m_numberFormat.format(geneStats.multiDipCombinations.get(dip)));
        rowNum += 1;
      }
    }
  }

  private void writeBgHeaders(Row headerRow, CellStyle headerStyle, int total, int geoTotal,
      SortedMap<String, Integer> regions) {
    if (doPivot()) {
      int cellCount = 2;
      if (total != geoTotal) {
        cellCount += 1;
        writeCell(headerRow, cellCount, "Samples with data: " + m_numberFormat.format(geoTotal), headerStyle);
      }
      for (String region : regions.keySet()) {
        cellCount += 1;
        writeCell(headerRow, cellCount, region + " Count (n=" +
            m_numberFormat.format(regions.get(region)) + ")", headerStyle);
        cellCount += 1;
        writeCell(headerRow, cellCount, region + " Frequency", headerStyle);
      }
    }
  }


  /**
   * Stringify the percentage so that we at least have 1 significant digit.
   */
  private String stringifyPercentage(int count, int total) {
    if (count == 0) {
      return "0%";
    }
    double percent = (double) count / total;
    // default is 3 decimal places
    if (percent >= 0.00001) {
      return m_percentFormat.format(percent);
    }
    String scientificNotation = Double.toString(percent);
    int idx = scientificNotation.indexOf("E-");
    // subtract 2 from exponent to account for the percentage shift (i.e. x100)
    int numZeros = Integer.parseInt(scientificNotation.substring(idx + 2)) - 2;
    m_percentFormat.setMaximumFractionDigits(numZeros);
    String p = m_percentFormat.format(percent);
    m_percentFormat.setMinimumFractionDigits(3);
    return p;
  }


  private void writeRow(Row row, int count, int total, int bgTotal, SortedMap<String, SortedMap<String, Integer>> bgMap,
      SortedMap<String, Integer> regions, String key) {

    writeCell(row, 0, key);
    writeCell(row, 1, m_numberFormat.format(count));
    writeCell(row, 2, stringifyPercentage(count, total));
    if (doPivot()) {
      int cellCount = 2;
      if (total != bgTotal) {
        cellCount += 1;
        writeCell(row, cellCount, "");
      }
      Map<String, Integer> pivotCounts = bgMap.get(key);
      for (String region : regions.keySet()) {
        int numForBg = 0;
        if (pivotCounts != null) {
          numForBg = pivotCounts.getOrDefault(region, 0);
        }
        cellCount += 1;
        writeCell(row, cellCount, m_numberFormat.format(numForBg));
        cellCount += 1;
        writeCell(row, cellCount, stringifyPercentage(numForBg, regions.get(region)));
      }
    }
  }


  private class GeneStats {
    private final DefinitionFile m_definitionFile;
    int numSamples = 0;
    int numUnknownSamples = 0;
    int numMultiDipSamples = 0;

    // calculated based on data
    int numSamplesComplete = 0;
    int numUnknownSamplesComplete = 0;
    int numMultiDipSamplesComplete = 0;

    int numSamplesMissingPositions = 0;
    int numSamplesUndocumentedVariants = 0;
    int numSamplesMissingAndUndocumented = 0;

    int numUnknownSamplesMissingPositions = 0;
    int numUnknownSamplesUndocumentedVariants = 0;
    int numUnknownSamplesMissingAndUndocumented = 0;

    int numMultiDipSamplesMissingPositions = 0;
    int numMultiDipSamplesUndocumentedVariants = 0;
    int numMultiDipSamplesMissingAndUndocumented = 0;


    int numAlleles = 0;
    SortedMap<String, Integer> diplotypeCounts = new TreeMap<>();
    SortedMap<String, Integer> phenotypeCounts = new TreeMap<>();
    SortedMap<String, Integer> alleleCounts = new TreeMap<>();

    SortedMap<String, Integer> missingPositions = new TreeMap<>();
    SortedMap<String, Integer> undocumentedVariants = new TreeMap<>();
    SortedMap<String, String> undocumentedRsids = new TreeMap<>();

    SortedMap<String, Integer> multiDipCombinations = new TreeMap<>();
    SortedMap<String, Integer> multiDipMissingPositions = new TreeMap<>();
    SortedMap<String, Integer> multiDipUndocumentedVariants = new TreeMap<>();

    SortedMap<String, Integer> unknownDipCombinations = new TreeMap<>();
    SortedMap<String, Integer> unknownMissingPositions = new TreeMap<>();
    SortedMap<String, Integer> unknownUndocumentedVariants = new TreeMap<>();

    int numBgSamples = 0;
    int numBgMultiDips = 0;
    int numBgAlleles = 0;
    SortedMap<String, SortedMap<String, Integer>> bgDiplotypeCounts = new TreeMap<>();
    SortedMap<String, SortedMap<String, Integer>> bgPhenotypeCounts = new TreeMap<>();
    SortedMap<String, SortedMap<String, Integer>> bgAlleleCounts = new TreeMap<>();
    SortedMap<String, Integer> diplotypeRegions = new TreeMap<>();
    SortedMap<String, Integer> phenotypeRegions = new TreeMap<>();
    SortedMap<String, Integer> alleleRegions = new TreeMap<>();


    private GeneStats(String gene) {
      m_definitionFile = m_env.getDefinitionReader().getDefinitionFile(gene);
    }

    private boolean hasPosition(@Nullable String value) {
      return !StringUtils.isBlank(value) && !"no".equals(value);
    }


    private static final Pattern sf_hgvsVariantPattern = Pattern.compile("g\\.([\\d_]+)([ACGT])+>([ACGT]+)$");

    /**
     * Find undocumented variants based on diplotype.
     */
    private boolean findUndocumentedVariants(String diplotypes) {
      // don't bother checking multiple diplotypes or haplotypes
      if (diplotypes.contains(" AND ")) {
        return false;
      }
      List<String> dips;
      if (diplotypes.contains(" OR ")) {
        dips = List.of(diplotypes.split(" OR "));
      } else {
        dips = List.of(diplotypes);
      }
      Set<String> undocumented = new HashSet<>();
      for (String dip : dips) {
        for (String hap : dip.split("/")) {
          List<String> alleles;
          if (hap.startsWith("[")) {
            hap = hap.substring(1, hap.length() - 1);
            alleles = List.of(hap.split(CombinationMatcher.COMBINATION_JOINER_REGEX));
          } else if (hap.equals(Haplotype.UNKNOWN)) {
            continue;
          } else {
            alleles = List.of(hap);
          }
          for (String a : alleles) {
            NamedAllele na = m_definitionFile.getNamedAllele(a);
            if (na == null) {
              Matcher m = sf_hgvsVariantPattern.matcher(a);
              if (m.matches()) {
                long pos = Long.parseLong(m.group(1));
                String ref = m.group(2);
                String alt = m.group(3);
                VariantLocus vl = m_definitionFile.getVariantForPosition(pos);
                if (ref.equals(vl.getRef())) {
                  if (!vl.getAlts().contains(alt)) {
                    undocumented.add(a);
                    if (vl.getRsid() != null) {
                      undocumentedRsids.put(a, vl.getRsid());
                    }
                  }
                } else {
                  undocumented.add(a);
                }
              } else {
                undocumented.add(a);
              }
            }
          }
        }
      }
      if (!undocumented.isEmpty()) {
        boolean unknown = diplotypes.contains("Unknown") || flipToUnknown(diplotypes, true);
        boolean multiDips = dips.size() > 1;
        for (String v : undocumented) {
          undocumentedVariants.put(v, undocumentedVariants.getOrDefault(v, 0) + 1);
          if (unknown) {
            unknownUndocumentedVariants.put(v, unknownUndocumentedVariants.getOrDefault(v, 0) + 1);
          } else if (multiDips) {
            multiDipUndocumentedVariants.put(v, multiDipUndocumentedVariants.getOrDefault(v, 0) + 1);
          }
        }
        return true;
      }
      return false;
    }


    private List<String> parsePositions(String value) {
      if ("yes".equals(value)) {
        return Collections.emptyList();
      }
      List<String> variants = new ArrayList<>();
      sf_commaSplitter.split(value).forEach(pos -> {
        VariantLocus vl = m_definitionFile.getVariantForPosition(Long.parseLong(pos));
        if (vl.getRsid() != null) {
          variants.add(vl.getRsid());
        } else {
          variants.add(vl.getVcfChrPosition());
        }
      });
      return variants;
    }

    private boolean flipToUnknown(String dip, boolean hasUndocumented) {
      if (m_treatCombinationsAsUnknown && !dip.contains("Unknown")) {
        if (hasUndocumented || dip.contains("[")) {
          return true;
        }
        if (!dip.contains(" AND ")) {
          // look for documented partials (AND calls only lists haplotypes)
          String[] dips = dip.split(" OR ");
          for (String d : dips) {
            for (String allele : d.split("/")) {
              if (m_definitionFile.getNamedAllele(allele) == null) {
                return true;
              }
            }
          }
        }
      }
      return false;
    }

    private void add(String dip, @Nullable String pheno, @Nullable String bgGroup, @Nullable String missingPositions) {
      numSamples += 1;

      boolean doAlleles = !dip.equals(CallsOnlyFormat.NO_CALL_TAG);
      boolean multiDip = dip.contains(" OR ");

      boolean hasUndocumented = findUndocumentedVariants(dip);
      boolean hasMissingPositions = false;
      if (hasPosition(missingPositions)) {
        hasMissingPositions = true;
        List<String> rsids = parsePositions(missingPositions);
        rsids.forEach(pos -> this.missingPositions.put(pos, this.missingPositions.getOrDefault(pos, 0) + 1));

        if (dip.contains("Unknown") || flipToUnknown(dip, hasUndocumented)) {
          rsids.forEach(pos -> this.unknownMissingPositions.put(pos,
              this.unknownMissingPositions.getOrDefault(pos, 0) + 1));
        } else if (multiDip) {
          rsids.forEach(pos -> this.multiDipMissingPositions.put(pos,
              this.multiDipMissingPositions.getOrDefault(pos, 0) + 1));
        }
      }

      if (flipToUnknown(dip, hasUndocumented)) {
        unknownDipCombinations.put(dip, unknownDipCombinations.getOrDefault(dip, 0) + 1);
        dip = "Unknown/Unknown";
        pheno = null;
        multiDip = false;
      }

      boolean isUnknown = dip.contains("Unknown");
      if (isUnknown) {
        numUnknownSamples += 1;
      } else if (multiDip) {
        numMultiDipSamples += 1;
        multiDipCombinations.put(dip, multiDipCombinations.getOrDefault(dip, 0) + 1);
        dip = "multiple calls";
        doAlleles = false;
        if (bgGroup != null) {
          numBgMultiDips += 1;
        }
      }

      if (hasUndocumented) {
        numSamplesUndocumentedVariants += 1;
        if (isUnknown) {
          numUnknownSamplesUndocumentedVariants += 1;
        } else if (multiDip) {
          numMultiDipSamplesUndocumentedVariants += 1;
        }
      }
      if (hasMissingPositions) {
        numSamplesMissingPositions += 1;
        if (isUnknown) {
          numUnknownSamplesMissingPositions += 1;
        } else if (multiDip) {
          numMultiDipSamplesMissingPositions += 1;
        }
      }
      if (hasUndocumented && hasMissingPositions) {
        numSamplesMissingAndUndocumented += 1;
        if (isUnknown) {
          numUnknownSamplesMissingAndUndocumented += 1;
        } else if (multiDip) {
          numMultiDipSamplesMissingAndUndocumented += 1;
        }
      } else if (!hasUndocumented && !hasMissingPositions) {
        numSamplesComplete += 1;
        if (isUnknown) {
          numUnknownSamplesComplete += 1;
        } else if (multiDip) {
          numMultiDipSamplesComplete += 1;
        }
      }


      diplotypeCounts.put(dip, diplotypeCounts.getOrDefault(dip, 0) + 1);
      if (pheno != null && !pheno.isEmpty()) {
        if (multiDip) {
          pheno = "Combination";
        }
      } else {
        pheno = "No phenotype";
      }
      phenotypeCounts.put(pheno, phenotypeCounts.getOrDefault(pheno, 0) + 1);

      if (bgGroup != null) {
        numBgSamples += 1;
        diplotypeRegions.put(bgGroup, diplotypeRegions.getOrDefault(bgGroup, 0) + 1);
        phenotypeRegions.put(bgGroup, phenotypeRegions.getOrDefault(bgGroup, 0) + 1);

        SortedMap<String, Integer> dipGeoMap = bgDiplotypeCounts.computeIfAbsent(dip, k -> new TreeMap<>());
        dipGeoMap.put(bgGroup, dipGeoMap.getOrDefault(bgGroup, 0) + 1);
        SortedMap<String, Integer> phenoGeoMap = bgPhenotypeCounts.computeIfAbsent(pheno, k -> new TreeMap<>());
        phenoGeoMap.put(bgGroup, phenoGeoMap.getOrDefault(bgGroup, 0) + 1);
      }

      if (doAlleles) {
        List<String> alleles = new ArrayList<>();

        if (dip.contains(" AND ")) {
          for (String part : dip.split(" AND ")) {
            Collections.addAll(alleles, part.split("/"));
          }

        } else {
          Collections.addAll(alleles, dip.split("/"));
          if (alleles.size() > 2) {
            throw new IllegalStateException("More than 2 alleles: " + dip);
          }
        }

        numAlleles += alleles.size();
        for (String allele : alleles) {
          alleleCounts.put(allele, alleleCounts.getOrDefault(allele, 0) + 1);
        }

        if (bgGroup != null) {
          numBgAlleles += alleles.size();
          alleleRegions.put(bgGroup, alleleRegions.getOrDefault(bgGroup, 0) + alleles.size());
          for (String allele : alleles) {
            SortedMap<String, Integer> geoMap = bgAlleleCounts.computeIfAbsent(allele, k -> new TreeMap<>());
            geoMap.put(bgGroup, geoMap.getOrDefault(bgGroup, 0) + 1);
          }
        }
      }
    }
  }
}
