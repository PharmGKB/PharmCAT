package org.pharmgkb.pharmcat.phenotype;

import java.util.Map;
import java.util.TreeMap;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.TextConstants;


/**
 * Utility methods to help deal with Phenotypes
 */
public class PhenotypeUtils {

  private static final Map<String,String> ABBREVIATIONS = new TreeMap<>();
  static {
    ABBREVIATIONS.put("pm", "Poor Metabolizer");
    ABBREVIATIONS.put("im", "Intermediate Metabolizer");
    ABBREVIATIONS.put("nm", "Normal Metabolizer");
    ABBREVIATIONS.put("em", "Normal Metabolizer"); // not a typo, extensives are now normals
    ABBREVIATIONS.put("um", "Ultrarapid Metabolizer");
  }

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
    String trimmed = StringUtils.strip(value.replaceAll("\\s\\s+", " "));
    // use a "z" in metablizer
    trimmed = trimmed.replaceAll("[Mm][Ee][Tt][Aa][Bb][Oo][Ll][Ii][ZzSs][Ee][Rr][Ss]?", "Metabolizer");
    String lowered = trimmed.toLowerCase();

    // if this is a known abbreviation, expand it
    if (ABBREVIATIONS.containsKey(lowered)) {
      return ABBREVIATIONS.get(lowered);
    }

    if (lowered.contains(TextConstants.INDETERMINATE.toLowerCase())) {
      return TextConstants.INDETERMINATE;
    }

    // if any common pheno name is in the string, use that
    if (lowered.contains("extensive")) {
      lowered = lowered.replaceAll("extensive", "normal"); // rename extensives to normals
    }
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
