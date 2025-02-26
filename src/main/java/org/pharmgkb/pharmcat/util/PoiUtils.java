package org.pharmgkb.pharmcat.util;

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
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFColor;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * This class contains utilities for working with POI.
 *
 * @author Mark Woon
 */
public class PoiUtils {
  public static String POI_GREY = "cccccc";
  public static String POI_GREEN = "70eb8d";
  public static String POI_ORANGE = "ffc870";
  public static String POI_RED = "fe938c";
  public static String POI_BLUE ="add8e6";
  public static String POI_YELLOW = "e6d8ad";


  public static Cell writeCell(Row row, int colNum, @Nullable String value) {
    return writeCell(row, colNum, value, null);
  }

  public static Cell writeCell(Row row, int colNum, int value) {
    Cell cell = row.createCell(colNum, CellType.NUMERIC);
    cell.setCellValue(value);
    return cell;
  }

  public static Cell writeCell(Row row, int colNum, @Nullable String value, @Nullable CellStyle cellStyle) {
    Cell cell = row.createCell(colNum, CellType.STRING);
    if (StringUtils.stripToNull(value) != null) {
      cell.setCellValue(value);
    }
    if (cellStyle != null) {
      cell.setCellStyle(cellStyle);
    }
    return cell;
  }


  public static CellStyle createCellStyle(Workbook workbook, String rgbColor, boolean boldFont) {

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

  public static CellStyle getHeaderCellStyle(Workbook wb) {
    return createCellStyle(wb, POI_GREY, true);
  }
}
