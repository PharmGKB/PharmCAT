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

  private static final String sf_drugNameTemplate = "<div class=\"drugName\"><div class=\"%s\">%s <a href=\"#%s\">%s</a></div></div>";
  private static final String sf_variantAlleleTemplate = "<td class=\"%s\">%s%s</td>";

  private static final String sf_iconHighlight = "<i class=\"fa fa-square highlight\" aria-hidden=\"true\" title=\"Consult guideline\"></i>";
  private static final String sf_iconPossibly = "<i class=\"fa fa-exclamation-triangle rxPossibly\" aria-hidden=\"true\" title=\"Possible Rx change\"></i>";
  private static final String sf_iconChange = "<i class=\"fa fa-times rxChange\" aria-hidden=\"true\" title=\"Rx change\"></i>";
  private static final String sf_iconNormal = "<i class=\"fa fa-circle normal\" aria-hidden=\"true\" title=\"No Rx change\"></i>";

  public static String drug(DrugLink drug) {
    String cn;
    String icon;

    cn = " " + drug.getRxClass();

    switch (drug.getRxClass()) {
      case "rxChange":
        icon = sf_iconChange;
        break;
      case "rxPossibly":
        icon = sf_iconPossibly;
        break;
      default:
        icon = sf_iconNormal;
        break;
    }
    if (drug.isHighlighted()) {
      cn = " highlightDrug";
      icon = sf_iconHighlight;
    }

    return String.format(sf_drugNameTemplate, cn, icon, drug.getGuidelineId(), drug.getName());
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
