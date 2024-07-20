package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.xml.bind.DatatypeConverter;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.compress.utils.FileNameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.poi.ss.usermodel.BorderStyle;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.CellType;
import org.apache.poi.ss.usermodel.FillPatternType;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import org.apache.poi.xssf.usermodel.XSSFColor;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.common.util.CliHelper;
import org.pharmgkb.pharmcat.ReportableException;
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
 * @author Mark Woon
 */
public class Subsetter {
  private final Multimap<String, String> m_allowList;
  private final SortedSet<GeneData> m_geneData = new TreeSet<>();
  private enum HeaderType { DESC, RSID, CHR_POS, SUBSECTION }


  public Subsetter(Multimap<String, String> allowList) {
    m_allowList = allowList;
  }


  /**
   * Subsets existing allele definition based on allowlist.
   */
  public void subset() throws IOException {

    DefinitionReader definitionReader = DefinitionReader.defaultReader();

    // make sure each gene has reference because we require it
    for (String gene : m_allowList.keySet()) {
      m_allowList.put(gene, definitionReader.getReferenceHaplotype(gene).getName());
    }

    for (String gene : m_allowList.keySet()) {
      DefinitionFile origDefinition = definitionReader.getDefinitionFile(gene);
      GeneData geneData = new GeneData(origDefinition);

      for (String haplotypeName : m_allowList.get(gene)) {
        NamedAllele hap = origDefinition.getNamedAllele(haplotypeName);
        if (hap == null) {
          System.out.println("WARNING: " + gene + " does not have an allele named \"" + haplotypeName + "\" - IGNORING!");
          continue;
        }
        geneData.haplotypes.add(hap);
        if (!hap.isReference()) {
          if (gene.equals("CYP2C19") && haplotypeName.equals("*1")) {
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
            .filter(na -> !na.isReference() && !(gene.equals("CYP2C19") && na.getName().equals("*1")))
            .allMatch(na -> {
          String allele = na.getCpicAllele(vl);
          return allele == null || Iupac.isWobble(allele);
        })) {
          wobblePositions.add(vl);
        }
      }
      if (!wobblePositions.isEmpty()) {
        System.out.println(gene + " has " + wobblePositions.size() + " wobble-only position" +
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
      m_geneData.add(geneData);
    }
  }


  /**
   * Exports allowed named alleles and positions for validation.
   */
  public void exportData(Path file) throws IOException {
    if (!Files.exists(file.getParent())) {
      Files.createDirectories(file.getParent());
    }

    PhenotypeMap phenotypeMap = new PhenotypeMap();
    try (XSSFWorkbook workbook = new XSSFWorkbook()) {
      List<CellStyle> cellStyles = new ArrayList<>();
      cellStyles.add(createCellStyle(workbook, "cccccc", true));
      cellStyles.add(createCellStyle(workbook, "70eb8d", true));
      cellStyles.add(createCellStyle(workbook, "70eb8d", false));
      cellStyles.add(createCellStyle(workbook, "ffc870", true));
      cellStyles.add(createCellStyle(workbook, "ffc870", false));
      cellStyles.add(createCellStyle(workbook, "fe938c", true));
      cellStyles.add(createCellStyle(workbook, "fe938c", false));

      for (GeneData gd : m_geneData) {
        Sheet sheet = workbook.createSheet(gd.gene);

        // position description headers
        Row row = sheet.createRow(0);
        writePositionHeaders(row, gd, HeaderType.DESC, cellStyles);

        // rsid headers
        row = sheet.createRow(1);
        writePositionHeaders(row, gd, HeaderType.RSID, cellStyles);

        // chr:pos headers
        row = sheet.createRow(2);
        writeCell(row, 0, "Allele", cellStyles.get(0));
        writeCell(row, 1, "CPIC", cellStyles.get(0));
        writeCell(row, 2, "DPWG", cellStyles.get(0));
        writePositionHeaders(row, gd, HeaderType.CHR_POS, cellStyles);


        // main section
        row = sheet.createRow(3);
        writeCell(row, 0, "CALLABLE ALLELES", cellStyles.get(1));
        writeCell(row, 1, null, cellStyles.get(1));
        writeCell(row, 2, null, cellStyles.get(1));
        updateSectionRowStyle(row, cellStyles, gd);
        int rowNum =  writeSection(sheet, 3, phenotypeMap, gd, gd.haplotypes);

        if (!gd.subsetHaplotypes.isEmpty()) {
          rowNum += 2;
          row = sheet.createRow(rowNum);
          writeCell(row, 0, "SUBSET ALLELES", cellStyles.get(1));
          writeCell(row, 1, "not called (alleles/nucleotides not included)", cellStyles.get(2));
          writeCell(row, 2, "but all positions are called", cellStyles.get(2));
          updateSectionRowStyle(row, cellStyles, gd);
          rowNum =  writeSection(sheet, rowNum, phenotypeMap, gd, gd.subsetHaplotypes);
        }

        if (!gd.overlapHaplotypes.isEmpty()) {
          rowNum += 2;
          row = sheet.createRow(rowNum);
          writeCell(row, 0, "OVERLAP ALLELES", cellStyles.get(3));
          writeCell(row, 1, "not called", cellStyles.get(4));
          writeCell(row, 2, "because not all positions are included", cellStyles.get(4));
          updateSectionRowStyle(row, cellStyles, gd);
          rowNum =  writeSection(sheet, rowNum, phenotypeMap, gd, gd.overlapHaplotypes);
        }

        if (!gd.unusedHaplotypes.isEmpty()) {
          rowNum += 2;
          row = sheet.createRow(rowNum);
          writeCell(row, 0, "UNUSED ALLELES", cellStyles.get(5));
          writeCell(row, 1, "not called", cellStyles.get(6));
          writeCell(row, 2, "because no positions are included", cellStyles.get(6));
          updateSectionRowStyle(row, cellStyles, gd);
          writeSection(sheet, rowNum, phenotypeMap, gd, gd.unusedHaplotypes);
        }

        int numColumns = 8 + gd.variants.size() + gd.missedVariants.size() + gd.unusedVariants.size();
        for (int x = 0; x <= numColumns; x += 1) {
          sheet.autoSizeColumn(x);
        }
      }

      System.out.println("Writing to " + file);
      try (OutputStream fos = Files.newOutputStream(file)) {
        workbook.write(fos);
      }
    }
  }

  private CellStyle createCellStyle(Workbook workbook, String rgbColor, boolean boldFont) {

    CellStyle cellStyle = workbook.createCellStyle();
    cellStyle.setFillForegroundColor(IndexedColors.LEMON_CHIFFON.getIndex());
    cellStyle.setFillForegroundColor(new XSSFColor(DatatypeConverter.parseHexBinary(rgbColor)));
    cellStyle.setFillPattern(FillPatternType.SOLID_FOREGROUND);
    cellStyle.setBottomBorderColor(IndexedColors.GREY_25_PERCENT.getIndex());
    cellStyle.setRightBorderColor(IndexedColors.GREY_25_PERCENT.getIndex());
    cellStyle.setTopBorderColor(IndexedColors.GREY_25_PERCENT.getIndex());
    cellStyle.setBorderBottom(BorderStyle.THIN);
    cellStyle.setBorderRight(BorderStyle.THIN);
    cellStyle.setBorderTop(BorderStyle.THIN);
    if (boldFont) {
      Font font = workbook.createFont();
      font.setBold(true);
      cellStyle.setFont(font);
    }
    return cellStyle;
  }

  private int writeSection(Sheet sheet, int rowNum, PhenotypeMap phenotypeMap, GeneData gd,
      SortedSet<NamedAllele> haplotypes) {

    for (NamedAllele hap : haplotypes) {
      rowNum += 1;
      Row row = sheet.createRow(rowNum);
      writeHaplotype(row, phenotypeMap, gd, hap);
    }
    return rowNum;
  }

  private Cell writeCell(Row row, int colNum, @Nullable String value) {
    return writeCell(row, colNum, value, null);
  }

  private Cell writeCell(Row row, int colNum, @Nullable String value, @Nullable CellStyle cellStyle) {
    Cell cell = row.createCell(colNum, CellType.STRING);
    if (StringUtils.stripToNull(value) != null) {
      cell.setCellValue(value);
    }
    if (cellStyle != null) {
      cell.setCellStyle(cellStyle);
    }
    return cell;
  }

  private void writePositionHeaders(Row row, GeneData gd, HeaderType headerType, List<CellStyle> cellStyles) {
    int colNum = writePositionHeader(row, 2, gd.variants, headerType, cellStyles.get(1));
    if (headerType == HeaderType.DESC) {
      writeCell(row, 3, "Included Positions")
          .setCellStyle(cellStyles.get(1));
    }
    if (!gd.missedVariants.isEmpty()) {
      int startCol = colNum;
      colNum = writePositionHeader(row, colNum + 2, gd.missedVariants, headerType, cellStyles.get(3));
      if (headerType == HeaderType.DESC) {
        writeCell(row, startCol + 3, "Missing Overlap Positions")
            .setCellStyle(cellStyles.get(3));
      }
    }
    if (!gd.unusedVariants.isEmpty()) {
      writePositionHeader(row, colNum + 2, gd.unusedVariants, headerType, cellStyles.get(5));
      if (headerType == HeaderType.DESC) {
        writeCell(row, colNum + 3, "Unused Positions")
            .setCellStyle(cellStyles.get(5));
      }
    }
  }

  private int writePositionHeader(Row row, int colNum, SortedSet<VariantLocus> variants, HeaderType headerType,
      CellStyle cellStyle) {
    for (VariantLocus vl : variants) {
      colNum += 1;
      if (headerType == HeaderType.RSID) {
        writeCell(row, colNum, vl.getRsid() != null ? vl.getRsid() : null)
            .setCellStyle(cellStyle);
      } else if (headerType == HeaderType.CHR_POS) {
        writeCell(row, colNum, vl.getVcfChrPosition())
            .setCellStyle(cellStyle);
      } else {
        writeCell(row, colNum, null)
            .setCellStyle(cellStyle);
      }
    }
    return colNum;
  }

  private void writeHaplotype(Row row, PhenotypeMap phenotypeMap, GeneData gd, NamedAllele hap) {
    writeCell(row, 0, hap.getName());
    writeFunction(row, phenotypeMap, gd.gene, hap, DataSource.CPIC);
    writeFunction(row, phenotypeMap, gd.gene, hap, DataSource.DPWG);

    int rowCol = writeVariants(row, hap, gd.variants, 2) + 2;
    if (!gd.missedVariants.isEmpty()) {
      rowCol = writeVariants(row, hap, gd.missedVariants, rowCol) + 2;
    }
    if (!gd.unusedVariants.isEmpty()) {
      writeVariants(row, hap, gd.unusedVariants, rowCol);
    }
  }

  private void writeFunction(Row row, PhenotypeMap phenotypeMap, String gene, NamedAllele hap, DataSource src) {
    GenePhenotype gp = phenotypeMap.getPhenotype(gene, src);
    if (gp != null) {
      StringBuilder builder = new StringBuilder()
          .append(gp.getHaplotypeFunction(hap.getName()));
      String activity = gp.getHaplotypeActivity(hap.getName());
      if (activity != null) {
        builder.append(" (")
            .append(activity)
            .append(")");
      }
      int col = switch (src) {
        case CPIC -> 1;
        case DPWG -> 2;
        default -> throw new UnsupportedOperationException("Cannot handle " + src);
      };
      writeCell(row, col, builder.toString());
    }
  }

  private int writeVariants(Row row, NamedAllele hap, SortedSet<VariantLocus> variants, int startCol) {
    for (VariantLocus vl : variants) {
      startCol += 1;
      String allele = hap.getCpicAllele(vl);
      if (allele != null) {
        writeCell(row, startCol, allele);
      }
    }
    return startCol;
  }

  private void updateSectionRowStyle(Row row, List<CellStyle> cellStyles, GeneData gd) {
    int endOfVariants = 3 + gd.variants.size();
    int startOfMissedVariants = endOfVariants + 1;
    int endOfMissedVariants = startOfMissedVariants + gd.missedVariants.size() + 1;
    if (gd.missedVariants.isEmpty()) {
      startOfMissedVariants = endOfVariants;
      endOfMissedVariants = endOfVariants;
    }
    int startOfUnusedVariants = endOfMissedVariants + 1;
    int endOfUnusedVariants = startOfUnusedVariants + gd.unusedVariants.size() + 1;
    if (gd.unusedVariants.isEmpty()) {
      startOfUnusedVariants = endOfMissedVariants;
      endOfUnusedVariants = endOfMissedVariants;
    }

    for (int x = 3; x < endOfUnusedVariants; x += 1) {
      Cell cell = row.getCell(x);
      if (cell == null) {
        cell = row.createCell(x, CellType.STRING);
      }
      if (x < endOfVariants) {
        cell.setCellStyle(cellStyles.get(2));
      } else if (x > startOfMissedVariants && x < endOfMissedVariants) {
        cell.setCellStyle(cellStyles.get(4));
      } else if (x > startOfUnusedVariants) {
        cell.setCellStyle(cellStyles.get(6));
      }
    }
  }


  /**
   * Exports updated definition files.
   */
  public void exportDefinitionFiles(Path dir) throws IOException {
    if (!Files.exists(dir)) {
      Files.createDirectories(dir);
    }

    System.out.println("Saving allele definitions in " + dir);

    DefinitionReader definitionReader = DefinitionReader.defaultReader();
    DataSerializer dataSerializer = new DataSerializer();
    Set<DefinitionExemption> exemptions = new HashSet<>();
    try (VcfHelper vcfHelper = new VcfHelper()) {
      for (GeneData gd : m_geneData) {
        gd.definitionFile.getNamedAlleles().clear();
        gd.definitionFile.getNamedAlleles().addAll(gd.haplotypes);

        SortedSet<VariantLocus> ignoredPositions = Arrays.stream(gd.definitionFile.getVariants())
            .filter(vl -> !gd.variants.contains(vl))
            .collect(Collectors.toCollection(TreeSet::new));
        InternalWrapper.removeIgnoredPositions(gd.definitionFile, ignoredPositions, false);

        InternalWrapper.doVcfTranslation(gd.definitionFile, vcfHelper);

        // export
        Path jsonFile = dir.resolve(gd.gene + "_translation.json");
        dataSerializer.serializeToJson(gd.definitionFile, jsonFile);
        System.out.println("Wrote " + jsonFile);

        DefinitionExemption exemption = definitionReader.getExemption(gd.gene);
        if (exemption != null) {
          exemptions.add(exemption);
        }
      }
    }

    Path exemptionsFile = dir.resolve(DataManager.EXEMPTIONS_JSON_FILE_NAME);;
    dataSerializer.serializeToJson(exemptions, exemptionsFile);

    DataManager.exportVcfData(dir);
  }



  public static void main(String[] args) {

    try {
      CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
          .addVersion("PharmCAT Definition Subsetter" + CliUtils.getVersion())
          .addOption("a", "alleles", "Allele whitelist file (tsv or xlsx)", false, "file")
          .addOption("p", "positions", "Positions whitelist file", false, "file")
          .addOption("o", "output-dir", "File containing a list of sample, one per line", false, "directory");
      if (!cliHelper.parse(args)) {
        return;
      }


      Path file = null;
      Subsetter subsetter = null;
      if (cliHelper.hasOption("p") && cliHelper.hasOption("a")) {
        System.out.println("-p and -a are mutually exclusive.");
        return;
      }

      if (cliHelper.hasOption("p")) {
        file = cliHelper.getValidFile("p", true);
        if (!file.toString().endsWith(".bim")) {
          throw new ReportableException("Unsupported file type for positions whitelist: " + file);
        }
        Multimap<String, String> allowList = parsePositionsBim(file);
        subsetter = new Subsetter(allowList);

      } else if (cliHelper.hasOption("a")) {
        file = cliHelper.getValidFile("a", true);
        Multimap<String, String> allowList;
        if (file.toString().endsWith(".tsv")) {
          allowList = parseAlleleAllowListTsv(file);
        } else if (file.toString().endsWith(".xlsx")) {
          allowList = parseAlleleAllowListXls(file);
        } else {
          throw new ReportableException("Unsupported file type for alleles whitelist: " + file);
        }
        subsetter = new Subsetter(allowList);
      }

      if (subsetter == null) {
        System.out.println("Nothing to do...");
        return;
      }

      subsetter.subset();

      Path dir;
      if (cliHelper.hasOption("o")) {
        dir = cliHelper.getValidDirectory("o", true);
      } else {
        dir = file.getParent();
      }

      String basename = FileNameUtils.getBaseName(file);
      subsetter.exportData(dir.resolve(basename + "-subset_details.xlsx"));
      subsetter.exportDefinitionFiles(dir.resolve(basename + "-subset_definitions"));

      System.out.println("Done.");

    } catch (CliHelper.InvalidPathException | ReportableException ex) {
      System.out.println(ex.getMessage());
    } catch (Exception e) {
      //noinspection CallToPrintStackTrace
      e.printStackTrace();
    }
  }


  private static Multimap<String, String> parsePositionsBim(Path file) throws IOException {

    DefinitionReader definitionReader = DefinitionReader.defaultReader();
    Set<String> rsidsOfInterest = definitionReader.getLocationsOfInterest().values().stream()
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
    Multimap<String, String> allowList = HashMultimap.create();
    for (String gene : definitionReader.getGenes()) {
      VariantLocus[] positions = definitionReader.getPositions(gene);

      for (NamedAllele na : definitionReader.getHaplotypes(gene)) {
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
          allowList.put(gene, na.getName());
        }
      }
    }

    return allowList;
  }


  private static Multimap<String, String> parseAlleleAllowListTsv(Path tsvFile) throws IOException {
    Multimap<String, String> allowList = HashMultimap.create();
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
        allowList.put(gene, allele);
        line = reader.readLine();
      }
    }
    return allowList;
  }

  private static Multimap<String, String> parseAlleleAllowListXls(Path xlsxFile) throws IOException {
    Multimap<String, String> allowList = HashMultimap.create();

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
        String allele = alleleCell.getStringCellValue();
        if (allele == null) {
          continue;
        }
        if (allele.startsWith("\"") && allele.endsWith("\"")) {
          allele = allele.substring(1, allele.length() - 1);
        }
        allowList.put(gene, allele);
      }
    }
    return allowList;
  }


  private static final class GeneData implements Comparable<GeneData> {
    private final DefinitionFile definitionFile;
    private final String gene;
    private final SortedSet<NamedAllele> haplotypes = new TreeSet<>();
    private final SortedSet<NamedAllele> subsetHaplotypes = new TreeSet<>();
    private final SortedSet<NamedAllele> overlapHaplotypes = new TreeSet<>();
    private final SortedSet<NamedAllele> unusedHaplotypes = new TreeSet<>();

    private final SortedSet<VariantLocus> variants = new TreeSet<>();
    private final SortedSet<VariantLocus> missedVariants = new TreeSet<>();
    private final SortedSet<VariantLocus> unusedVariants = new TreeSet<>();


    private GeneData(DefinitionFile definitionFile) {
      this.definitionFile = definitionFile;
      gene = definitionFile.getGeneSymbol();
    }

    private void calculateUnused() {
      Arrays.stream(definitionFile.getVariants())
          .filter(vl -> !variants.contains(vl) && !missedVariants.contains(vl))
          .forEach(unusedVariants::add);
      definitionFile.getNamedAlleles().stream()
          .filter(na -> !haplotypes.contains(na) && !subsetHaplotypes.contains(na) && !overlapHaplotypes.contains(na))
          .forEach(unusedHaplotypes::add);
    }

    @Override
    public int hashCode() {
      return Objects.hashCode(gene);
    }

    @Override
    public boolean equals(Object o) {
      if (o == this) {
        return true;
      }
      if (!(o instanceof GeneData gd)) {
        return false;
      }
      return Objects.equals(gene, gd.gene);
    }

    @Override
    public int compareTo(GeneData o) {
      return gene.compareTo(o.gene);
    }
  }
}
