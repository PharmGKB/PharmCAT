package org.pharmgkb.pharmcat.reporter.model;

/**
 * Data source names.
 */
public enum DataSource {
  // XXX: code expects CPIC to be sorted first (e.g. in HtmlFormat.compile())
  // DO NOT CHANGE ORDER without dealing with this expectation
  CPIC("CPIC", "CPIC"),
  DPWG("DPWG", "DPWG"),
  PHARMGKB("PharmGKB", "PharmGKB"),
  FDA("FDA", "FDA"),
  UNKNOWN("Unknown", "Unknown");

  private final String displayName;
  private final String pharmgkbName;

  DataSource(String displayName, String pharmgkbName) {
    this.displayName = displayName;
    this.pharmgkbName = pharmgkbName;
  }

  public String getDisplayName() {
    return this.displayName;
  }

  public String getPharmgkbName() {
    return this.pharmgkbName;
  }
}
