/*
 ----- BEGIN LICENSE BLOCK -----
 This Source Code Form is subject to the terms of the Mozilla Public License, v.2.0.
 If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
 ----- END LICENSE BLOCK -----
 */
package org.pharmgkb.common.io;

import java.io.IOException;
import java.io.PrintWriter;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.security.GeneralSecurityException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import com.google.gdata.client.spreadsheet.SpreadsheetService;
import com.google.gdata.data.spreadsheet.CellEntry;
import com.google.gdata.data.spreadsheet.CellFeed;
import com.google.gdata.data.spreadsheet.SpreadsheetEntry;
import com.google.gdata.data.spreadsheet.SpreadsheetFeed;
import com.google.gdata.data.spreadsheet.WorksheetEntry;
import com.google.gdata.util.ServiceException;


/**
 * This is a helper class for working with the Google Sheets API.
 * <p>
 * This relies on the old GData API.
 * Docs at https://developers.google.com/google-apps/spreadsheets/.
 *
 * @author Mark Woon
 */
public class GoogleSheetsHelper implements AutoCloseable {
  public static final String SHEETS_SCOPE = "https://spreadsheets.google.com/feeds";
  /** Pattern to parse R1C1 notation for row/column numbers. */
  private static final Pattern sf_rcPattern = Pattern.compile(".*/R(\\d+)C(\\d+)$");
  private final URL m_spreadsheetUrl = new URL("https://spreadsheets.google.com/feeds/spreadsheets/private/full");
  private boolean m_isPersonalApiHelper;
  private GoogleApiHelper m_googleApiHelper;
  private SpreadsheetService m_sheetsService;


  public GoogleSheetsHelper(@Nonnull String userId, @Nonnull String privateKey, @Nonnull String service)
      throws IOException, GeneralSecurityException {

    m_googleApiHelper = new GoogleApiHelper(userId, privateKey, SHEETS_SCOPE);
    m_isPersonalApiHelper = true;
    m_sheetsService = new SpreadsheetService(service);
    m_sheetsService.setOAuth2Credentials(m_googleApiHelper.getCredential());
  }


  public GoogleSheetsHelper(@Nonnull GoogleApiHelper apiHelper, @Nonnull String service)
      throws IOException, GeneralSecurityException {

    m_googleApiHelper = apiHelper;
    m_sheetsService = new SpreadsheetService(service);
    m_sheetsService.setOAuth2Credentials(m_googleApiHelper.getCredential());
  }


  @Override
  public void close() {
    if (m_isPersonalApiHelper) {
      m_googleApiHelper.close();
    }
  }


  public SpreadsheetService getService() {
    return m_sheetsService;
  }


  /**
   * Gets a list of all sheets for user.
   */
  public List<SpreadsheetEntry> getSheets() throws IOException, ServiceException {
    return m_sheetsService.getFeed(m_spreadsheetUrl, SpreadsheetFeed.class).getEntries();
  }


  /**
   * Exports sheet with the specified ID to TSV.
   *
   * @param fileId a file Id from Google Drive
   */
  public void exportToTsv(String fileId, Path tsvFile) throws IOException, ServiceException {
    exportToTsv(new URL("https://spreadsheets.google.com/feeds/spreadsheets/" + fileId), tsvFile);
  }

  /**
   * Exports sheet with the specified ID to TSV.
   *
   * @param url URL from {@link SpreadsheetEntry#getId()} ()}
   */
  public void exportToTsv(URL url, Path tsvFile) throws IOException, ServiceException {
    // fetch sheet
    SpreadsheetEntry spreadsheet = m_sheetsService.getEntry(url, SpreadsheetEntry.class);
    // get first worksheet
    WorksheetEntry worksheet = spreadsheet.getDefaultWorksheet();

    // fetch the cell feed of the worksheet
    URL cellFeedUrl = worksheet.getCellFeedUrl();
    CellFeed cellFeed = m_sheetsService.getFeed(cellFeedUrl, CellFeed.class);

    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(tsvFile))) {
      int curRow = 1;
      int curCol = 1;
      for (CellEntry cell : cellFeed.getEntries()) {
        Matcher m = sf_rcPattern.matcher(cell.getId());
        if (!m.matches()) {
          throw new IOException("Unexpected cell id (" + cell.getId() + ").  Cannot determine row/column.");
        }
        int row = Integer.parseInt(m.group(1));
        int col = Integer.parseInt(m.group(2));

        while (curRow < row) {
          writer.println();
          curRow += 1;
          curCol = 1;
        }
        while (curCol < col) {
          writer.print("\t");
          curCol += 1;
        }
        writer.print(cell.getCell().getValue());
      }
      writer.println();
    }
  }
}
