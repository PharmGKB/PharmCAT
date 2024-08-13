package org.pharmgkb.pharmcat.subsetter;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
import javax.xml.bind.DatatypeConverter;
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
import org.apache.poi.xssf.usermodel.XSSFColor;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.phenotype.PhenotypeMap;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.reporter.model.DataSource;


/**
 * @author Mark Woon
 */
public class SummaryWriter {
  private static final int CELL_STYLE_REPORT_AS_REF = 7;
  private static final int CELL_STYLE_MODIFIED_FUNCTION = 8;
  private enum HeaderType { DESC, RSID, CHR_POS, SUBSECTION }

  private final PhenotypeMap m_phenotypeMap;
  private final SortedMap<String, GeneData> m_geneData;

  private final List<CellStyle> m_cellStyles = new ArrayList<>();


  public SummaryWriter(PhenotypeMap phenotypeMap, SortedMap<String, GeneData> geneData) {
    m_phenotypeMap = phenotypeMap;
    m_geneData = geneData;
  }


  /**
   * Exports allowed named alleles and positions for validation.
   */
  public void write(Path file) throws IOException {
    if (!Files.exists(file.getParent())) {
      Files.createDirectories(file.getParent());
    }

    try (XSSFWorkbook workbook = new XSSFWorkbook()) {
      m_cellStyles.clear();
      m_cellStyles.add(createCellStyle(workbook, "cccccc", true));
      m_cellStyles.add(createCellStyle(workbook, "70eb8d", true));
      m_cellStyles.add(createCellStyle(workbook, "70eb8d", false));
      m_cellStyles.add(createCellStyle(workbook, "ffc870", true));
      m_cellStyles.add(createCellStyle(workbook, "ffc870", false));
      m_cellStyles.add(createCellStyle(workbook, "fe938c", true));
      m_cellStyles.add(createCellStyle(workbook, "fe938c", false));
      // for report as reference
      m_cellStyles.add(createCellStyle(workbook, "add8e6", false));
      // for modified function
      m_cellStyles.add(createCellStyle(workbook, "e6d8ad", false));

      for (String gene : m_geneData.keySet()) {
        GeneData gd = m_geneData.get(gene);
        Sheet sheet = workbook.createSheet(gd.gene);

        // position description headers
        Row row = sheet.createRow(0);
        writePositionHeaders(row, gd, HeaderType.DESC);

        // rsid headers
        row = sheet.createRow(1);
        writePositionHeaders(row, gd, HeaderType.RSID);

        // chr:pos headers
        row = sheet.createRow(2);
        writeCell(row, 0, "Allele", m_cellStyles.get(0));
        writeCell(row, 1, "CPIC", m_cellStyles.get(0));
        writeCell(row, 2, "DPWG", m_cellStyles.get(0));
        writePositionHeaders(row, gd, HeaderType.CHR_POS);


        // main section
        row = sheet.createRow(3);
        writeCell(row, 0, "CALLABLE ALLELES", m_cellStyles.get(1));
        writeCell(row, 1, null, m_cellStyles.get(1));
        writeCell(row, 2, null, m_cellStyles.get(1));
        updateSectionRowStyle(row, gd);
        int rowNum =  writeSection(sheet, 3, gd, gd.haplotypes);

        if (!gd.subsetHaplotypes.isEmpty()) {
          rowNum += 2;
          row = sheet.createRow(rowNum);
          writeCell(row, 0, "SUBSET ALLELES", m_cellStyles.get(1));
          writeCell(row, 1, "not called (alleles/nucleotides not included)", m_cellStyles.get(2));
          writeCell(row, 2, "but all positions are called", m_cellStyles.get(2));
          updateSectionRowStyle(row, gd);
          rowNum =  writeSection(sheet, rowNum, gd, gd.subsetHaplotypes);
        }

        if (!gd.overlapHaplotypes.isEmpty()) {
          rowNum += 2;
          row = sheet.createRow(rowNum);
          writeCell(row, 0, "OVERLAP ALLELES", m_cellStyles.get(3));
          writeCell(row, 1, "not called", m_cellStyles.get(4));
          writeCell(row, 2, "because not all positions are included", m_cellStyles.get(4));
          updateSectionRowStyle(row, gd);
          rowNum =  writeSection(sheet, rowNum, gd, gd.overlapHaplotypes);
        }

        if (!gd.unusedHaplotypes.isEmpty()) {
          rowNum += 2;
          row = sheet.createRow(rowNum);
          writeCell(row, 0, "UNUSED ALLELES", m_cellStyles.get(5));
          writeCell(row, 1, "not called", m_cellStyles.get(6));
          writeCell(row, 2, "because no positions are included", m_cellStyles.get(6));
          updateSectionRowStyle(row, gd);
          writeSection(sheet, rowNum, gd, gd.unusedHaplotypes);
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

  private int writeSection(Sheet sheet, int rowNum, GeneData gd, SortedSet<NamedAllele> haplotypes) {

    for (NamedAllele hap : haplotypes) {
      rowNum += 1;
      Row row = sheet.createRow(rowNum);
      writeHaplotype(row, gd, hap);
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

  private void writePositionHeaders(Row row, GeneData gd, HeaderType headerType) {
    int colNum = writePositionHeader(row, 2, gd.variants, headerType, m_cellStyles.get(1));
    if (headerType == HeaderType.DESC) {
      writeCell(row, 3, "Included Positions")
          .setCellStyle(m_cellStyles.get(1));
    }
    if (!gd.missedVariants.isEmpty()) {
      int startCol = colNum;
      colNum = writePositionHeader(row, colNum + 2, gd.missedVariants, headerType, m_cellStyles.get(3));
      if (headerType == HeaderType.DESC) {
        writeCell(row, startCol + 3, "Missing Overlap Positions")
            .setCellStyle(m_cellStyles.get(3));
      }
    }
    if (!gd.unusedVariants.isEmpty()) {
      writePositionHeader(row, colNum + 2, gd.unusedVariants, headerType, m_cellStyles.get(5));
      if (headerType == HeaderType.DESC) {
        writeCell(row, colNum + 3, "Unused Positions")
            .setCellStyle(m_cellStyles.get(5));
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

  private void writeHaplotype(Row row, GeneData gd, NamedAllele hap) {
    CellStyle cellStyle = null;
    writeCell(row, 0, hap.getName(), cellStyle);
    writeFunction(row, gd.gene, hap, DataSource.CPIC);
    writeFunction(row, gd.gene, hap, DataSource.DPWG);

    int rowCol = writeVariants(row, hap, gd.variants, 2) + 2;
    if (!gd.missedVariants.isEmpty()) {
      rowCol = writeVariants(row, hap, gd.missedVariants, rowCol) + 2;
    }
    if (!gd.unusedVariants.isEmpty()) {
      writeVariants(row, hap, gd.unusedVariants, rowCol);
    }
  }

  private void writeFunction(Row row, String gene, NamedAllele hap, DataSource src) {
    GenePhenotype gp = m_phenotypeMap.getPhenotype(gene, src);
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
      CellStyle cellStyle = gp.isModified(hap.getName()) ? m_cellStyles.get(CELL_STYLE_MODIFIED_FUNCTION) : null;
      writeCell(row, col, builder.toString(), cellStyle);
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

  private void updateSectionRowStyle(Row row, GeneData gd) {
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
        cell.setCellStyle(m_cellStyles.get(2));
      } else if (x > startOfMissedVariants && x < endOfMissedVariants) {
        cell.setCellStyle(m_cellStyles.get(4));
      } else if (x > startOfUnusedVariants) {
        cell.setCellStyle(m_cellStyles.get(6));
      }
    }
  }
}
