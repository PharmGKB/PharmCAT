package org.pharmgkb.pharmcat.reporter.model;

import java.util.Arrays;
import java.util.List;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.reporter.model.pgkb.DosingGuideline;


/**
 * This enum represents the different sources of prescribing guidance
 */
public enum PrescribingGuidanceSource {
  CPIC_GUIDELINE("CPIC Guideline Annotation", "cpic-guideline", DataSource.CPIC, DataSource.CPIC, "Guideline Annotation"),
  DPWG_GUIDELINE("DPWG Guideline Annotation", "dpwg-guideline", DataSource.DPWG, DataSource.DPWG, "Guideline Annotation"),
  FDA_LABEL("FDA Label Annotation", "fda-label", DataSource.CPIC, DataSource.FDA, "Label Annotation"),
  FDA_ASSOC("FDA PGx Association", "fda-assoc", DataSource.CPIC, DataSource.FDA, "PGx Association");

  private final String displayName;
  private final String codeName;
  private final DataSource phenoSource;
  private final DataSource pgkbSource;
  private final String pgkbObjectType;

  PrescribingGuidanceSource(String displayName, String codeName, DataSource phenoSource, DataSource pgkbSource, String pgkbObjectType) {
    this.displayName = displayName;
    this.codeName = codeName;
    this.phenoSource = phenoSource;
    this.pgkbSource = pgkbSource;
    this.pgkbObjectType = pgkbObjectType;
  }

  /**
   * Lists all available prescribing guidance sources
   */
  public static List<PrescribingGuidanceSource> listValues() {
    return Arrays.asList(values());
  }

  /**
   * Gets the text used to label this source in public-facing documentation
   */
  public String getDisplayName() {
    return this.displayName;
  }

  /**
   * Gets the string to use in code or structured text for this source
   */
  public String getCodeName() {
    return this.codeName;
  }

  /**
   * Gets the {@link DataSource} this uses to get phenotypes to match to recommendations
   */
  public DataSource getPhenoSource() {
    return this.phenoSource;
  }

  /**
   * Gets the {@link DataSource} this has been assigned in ClinPGx.
   */
  public DataSource getPgkbSource() {
    return this.pgkbSource;
  }

  /**
   * Gets the object type this has been assigned in ClinPGx.
   */
  public String getPgkbObjectType() {
    return this.pgkbObjectType;
  }

  public boolean matches(DosingGuideline prescribingGuidanceDocument) {
    return
        // matches the source
        prescribingGuidanceDocument.getSource().equals(getPgkbSource().getPharmgkbName())
        // and the annotation type
        && prescribingGuidanceDocument.getObjCls().equals(getPgkbObjectType());
  }

  public static @Nullable PrescribingGuidanceSource typeFor(DosingGuideline prescribingGuidanceDocument) {
    for (PrescribingGuidanceSource source : values()) {
      if (source.matches(prescribingGuidanceDocument)) {
        return source;
      }
    }
    return null;
  }

  public String toString() {
    return getDisplayName();
  }
}
