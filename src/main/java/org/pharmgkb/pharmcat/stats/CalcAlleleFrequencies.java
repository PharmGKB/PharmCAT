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
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.reporter.format.CallsOnlyFormat;
import org.pharmgkb.pharmcat.util.CliUtils;
import org.pharmgkb.pharmcat.util.PoiUtils;

import static org.pharmgkb.pharmcat.util.PoiUtils.writeCell;


/**
 * This tool generates allele frequency stats from *.report.tsv files.
 *
 * @author Mark Woon
 */
public class CalcAlleleFrequencies {
  private final int m_pivotCol;
  private final Map<String, GeneStats> m_stats = new TreeMap<>();
  private final NumberFormat m_numberFormat = NumberFormat.getNumberInstance();
  private final NumberFormat m_percentFormat = NumberFormat.getPercentInstance();


  public static void main(@Nullable String[] args) throws Exception {
    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addOption("i", "input", "Calls-only TSV file or directory", true, "file")
          .addOption("o", "output-dir", "Output to (optional, default is input file directory)", false, "directory")
          .addOption("pc", "pivot-column", "Pivot column number; first column is 1", false, "number")
          ;
      if (!cliHelper.parse(args)) {
        CliUtils.failIfNotTest();
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

      CalcAlleleFrequencies calcAlleleFrequencies = new CalcAlleleFrequencies(pivotCol);
      if (Files.isDirectory(input)) {
        calcAlleleFrequencies.ingestDir(input);
      } else {
        calcAlleleFrequencies.ingestFile(input);
      }
      calcAlleleFrequencies.write(outDir);

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      CliUtils.failIfNotTest(ex.getMessage());
    } catch (Exception ex) {
      CliUtils.failIfNotTest(ex);
    }
  }


  private CalcAlleleFrequencies(int pivotCol) {
    m_pivotCol = pivotCol;
    m_percentFormat.setMinimumFractionDigits(3);
  }

  private boolean doPivot() {
    return m_pivotCol >= 0;
  }


  private void ingestDir(Path dir) throws IOException {
    try (DirectoryStream<Path> stream = Files.newDirectoryStream(dir)) {
      for (Path path : stream) {
        if (Files.isDirectory(path)) {
          ingestDir(path);
        } else if (Files.isRegularFile(path) && path.toString().endsWith(".report.tsv")) {
          ingestFile(path);
        }
      }
    }
  }

  private void ingestFile(Path file) throws IOException {

    int geneCol = 0;
    System.out.println("Reading " + file);
    try (BufferedReader reader = Files.newBufferedReader(file)) {
      String line = reader.readLine();
      int lineNum = 1;
      // skip PharmCAT version info
      if (line.startsWith("PharmCAT")) {
        line = reader.readLine();
        lineNum += 1;
      }
      // check header line
      if (line.startsWith("Sample ID")) {
        geneCol = 1;
      } else if (!line.startsWith("Gene\t")) {
        throw new IOException("Expecting headers but line " + lineNum + " is: " + line);
      }

      while ((line = reader.readLine()) != null) {
        lineNum += 1;
//        if (lineNum % 5000 == 0) {
//          System.out.println(m_numberFormat.format(lineNum));
//        }
        String[] fields = line.split("\t");
        String gene = fields[geneCol];
        String dip = fields[geneCol + 1];
        String bgData = null;
        if (doPivot()) {
          if (m_pivotCol < fields.length) {
            bgData = fields[m_pivotCol];
          }
        }
        try {
          GeneStats geneStats = m_stats.computeIfAbsent(gene, k -> new GeneStats());
          geneStats.add(dip, bgData);
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

      writeRow(row, geneStats.diplotypeCounts.get(key), geneStats.numSamples, geneStats.numBgSamples, geneStats,
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
      writeRow(row, geneStats.alleleCounts.get(key), geneStats.numAlleles, geneStats.numBgAlleles, geneStats,
          geneStats.bgAlleleCounts, geneStats.alleleRegions, key);
    }
    for (int x = 0; x <= 2; x += 1) {
      sheet.autoSizeColumn(x);
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


  private void writeRow(Row row, int count, int total, int bgTotal, GeneStats geneStats,
      SortedMap<String, SortedMap<String, Integer>> bgMap, SortedMap<String, Integer> regions, String key) {

    writeCell(row, 0, key);
    writeCell(row, 1, m_numberFormat.format(count));
    writeCell(row, 2, m_percentFormat.format(count / (double) total));
    if (doPivot()) {
      int cellCount = 2;
      if (total != bgTotal) {
        cellCount += 1;
        writeCell(row, cellCount, "");
      }
      Map<String, Integer> pivotCounts = bgMap.get(key);
      for (String region : regions.keySet()) {
        int numDipsForBg = 0;
        if (pivotCounts != null) {
          numDipsForBg = pivotCounts.getOrDefault(region, 0);
        }
        cellCount += 1;
        writeCell(row, cellCount, m_numberFormat.format(numDipsForBg));
        cellCount += 1;
        writeCell(row, cellCount, m_percentFormat.format(numDipsForBg / (double)regions.get(region)));
      }
    }
  }


  private static class GeneStats {
    int numSamples = 0;
    int numMultiDips = 0;
    int numAlleles = 0;
    SortedMap<String, Integer> diplotypeCounts = new TreeMap<>();
    SortedMap<String, Integer> alleleCounts = new TreeMap<>();

    int numBgSamples = 0;
    int numBgMultiDips = 0;
    int numBgAlleles = 0;
    SortedMap<String, SortedMap<String, Integer>> bgDiplotypeCounts = new TreeMap<>();
    SortedMap<String, SortedMap<String, Integer>> bgAlleleCounts = new TreeMap<>();
    SortedMap<String, Integer> diplotypeRegions = new TreeMap<>();
    SortedMap<String, Integer> alleleRegions = new TreeMap<>();


    private void add(String dip, @Nullable String bgGroup) {
      numSamples += 1;
      boolean doAlleles = !dip.equals(CallsOnlyFormat.NO_CALL_TAG);

      if (dip.contains(" OR ")) {
        numMultiDips += 1;
        dip = "multiple calls";
        doAlleles = false;
        if (bgGroup != null) {
          numBgMultiDips += 1;
        }
      }

      diplotypeCounts.put(dip, diplotypeCounts.getOrDefault(dip, 0) + 1);

      if (bgGroup != null) {
        numBgSamples += 1;
        diplotypeRegions.put(bgGroup, diplotypeRegions.getOrDefault(bgGroup, 0) + 1);

        SortedMap<String, Integer> geoMap = bgDiplotypeCounts.computeIfAbsent(dip, k -> new TreeMap<>());
        geoMap.put(bgGroup, geoMap.getOrDefault(bgGroup, 0) + 1);
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
