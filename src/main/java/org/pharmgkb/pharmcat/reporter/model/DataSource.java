package org.pharmgkb.pharmcat.reporter.model;

/**
 * Data source names
 */
public enum DataSource {
  CPIC("CPIC", "CPIC"),
  // NOTE: we want to indicate that the Dutch (DPWG) data comes through PharmGKB annotations so we add it to the display
  DPWG("PharmGKB-DPWG", "DPWG"),
  PHARMGKB("PharmGKB", "PharmGKB"),
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
