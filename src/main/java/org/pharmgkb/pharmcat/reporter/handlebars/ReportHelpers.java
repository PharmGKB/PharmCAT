package org.pharmgkb.pharmcat.reporter.handlebars;

import org.pharmgkb.pharmcat.reporter.model.DrugLink;


/**
 * Class to hold methods that can be used as helper methods in the Handlebars templating system
 *
 * @author Ryan Whaley
 */
public class ReportHelpers {

  private static final String sf_drugNameTemplate = "<div class=\"%s\"><a href=\"#%s\">%s</a></div>";

  public static String drug(DrugLink drug) {
    String cn = "drugName";

    if (drug.isHighlighted()) {
      cn += " highlightDrug";
    }

    if (drug.isRxChange()) {
      cn += " rxChange";
    }

    return String.format(sf_drugNameTemplate, cn, drug.getGuidelineId(), drug.getName());
  }
}
