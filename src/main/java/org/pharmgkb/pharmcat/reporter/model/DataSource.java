package org.pharmgkb.pharmcat.reporter.model;

/**
 * Data source names
 */
public enum DataSource {
  CPIC("CPIC"),
  // NOTE: we want to indicate that the Dutch (DPWG) data comes through PharmGKB annotations so we add it to the display
  DPWG("PharmGKB-DPWG"),
  UNKNOWN("Unknown");

  private final String displayName;

  DataSource(String displayName) {
    this.displayName = displayName;
  }

  public String getDisplayName() {
    return this.displayName;
  }
}
