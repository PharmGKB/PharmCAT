package org.pharmgkb.pharmcat.reporter.handlebars;

import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.reporter.model.DrugLink;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;


/**
 * Class to hold methods that can be used as helper methods in the Handlebars templating system
 *
 * @author Ryan Whaley
 */
public class ReportHelpers {

  private static final String sf_drugNameTemplate = "<div class=\"drugName\"><a href=\"#%s\">%s</a></div>";
  private static final String sf_variantAlleleTemplate = "<td class=\"%s\">%s%s</td>";

  public static String drug(DrugLink drug) {
    return String.format(sf_drugNameTemplate, drug.getGuidelineId(), drug.getName());
  }

  public static String variantAlleles(VariantReport variantReport) {
    String cellStyle = variantReport.isNonwildtype() ? "nonwild" : "";
    String mismatch = variantReport.isMismatch() ? "<div class=\"callMessage\">Mismatch: Called allele does not match allele definitions</div>" : "";
    if (variantReport.isMismatch()) {
      cellStyle = StringUtils.strip(cellStyle + " mismatch");
    }
    return String.format(sf_variantAlleleTemplate, cellStyle, variantReport.getCall(), mismatch);
  }
}
