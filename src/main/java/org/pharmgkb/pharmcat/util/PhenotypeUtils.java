package org.pharmgkb.pharmcat.util;

import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import org.apache.commons.lang3.StringUtils;


/**
 * Utility methods to help deal with Phenotypes
 */
public class PhenotypeUtils {

  private static final Map<String,String> ABBREVIATIONS = new TreeMap<>();
  static {
    ABBREVIATIONS.put("PM", "Poor Metabolizer");
    ABBREVIATIONS.put("IM", "Intermediate Metabolizer");
    ABBREVIATIONS.put("NM", "Normal Metabolizer");
    ABBREVIATIONS.put("EM", "Normal Metabolizer"); // not a typo, extensives are now normals
    ABBREVIATIONS.put("UM", "Ultrarapid Metabolizer");
  }
  private static final String INDETERMINATE = "Indeterminate";

  /**
   * Normalize some common phenotype name strings into ones that PharmCAT understands
   * @param value a phenotype value
   * @return a normalized version of the phenotype value
   */
  public static String normalize(String value) {
    if (StringUtils.isBlank(value)) {
      return null;
    }

    // get rid of redundant whitespace
    String trimmed = StringUtils.trim(value.replaceAll("\\s\\s+", " "));

    // if this is a known abbreviation, expand it
    if (ABBREVIATIONS.containsKey(trimmed.toUpperCase(Locale.US))) {
      return ABBREVIATIONS.get(trimmed.toUpperCase(Locale.US));
    }

    String lowered = trimmed.toLowerCase();
    if (lowered.contains("etaboliser")) {
      lowered = lowered.replaceAll("etaboliser", "etabolizer"); // use a "z" in metablizer
    }

    if (lowered.contains("extensive")) {
      lowered = lowered.replaceAll("extensive", "normal"); // rename extensives to normals
    }

    if (lowered.contains("indeterminate")) {
      return INDETERMINATE;
    }

    // if any common pheno name is in the string, use that
    for (String phenoName : ABBREVIATIONS.values()) {
      if (lowered.contains(phenoName.toLowerCase())) {
        if (lowered.startsWith("likely")) {
          return "Likely " + phenoName;
        }
        if (lowered.startsWith("possible")) {
          return "Possible " + phenoName;
        }
        return phenoName;
      }
    }

    // if we can't find any common changes, just use the trimmed string
    return trimmed;
  }
}
