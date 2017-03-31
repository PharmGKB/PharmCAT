package org.pharmgkb.pharmcat.reporter.handlebars;

import org.pharmgkb.pharmcat.reporter.model.DrugLink;


/**
 * Class to hold methods that can be used as helper methods in the Handlebars templating system
 *
 * @author Ryan Whaley
 */
public class ReportHelpers {

  private static final String sf_drugNameTemplate = "<div class=\"drugName\"><div class=\"%s\"><a href=\"#%s\">%s</a> %s</div></div>";

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

    return String.format(sf_drugNameTemplate, cn, drug.getGuidelineId(), drug.getName(), icon);
  }

  public static String iconNormal() {
    return sf_iconNormal;
  }

  public static String iconPossibly() {
    return sf_iconPossibly;
  }

  public static String iconChange() {
    return sf_iconChange;
  }

  public static String iconHighlight() {
    return sf_iconHighlight;
  }
}
