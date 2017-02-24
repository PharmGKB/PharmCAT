package org.pharmgkb.pharmcat.reporter.handlebars;

import java.util.List;
import com.github.jknack.handlebars.Options;
import com.google.common.collect.ImmutableList;


/**
 * Class to hold methods that can be used as helper methods in the Handlebars templating system
 *
 * @author Ryan Whaley
 */
public class ReportHelpers {

  private static final String sf_drugNameTemplate = "<div class=\"%s\">%s</div>";
  private static final List<String> sf_highlightDrugs = ImmutableList.of(
      "warfarin"
  );

  public static String drug(String name, Options options) {
    String cn = "drugName";

    if (sf_highlightDrugs.contains(name)) {
      cn += " highlightDrug";
    }

    boolean rxChange = options.param(0);
    if (rxChange) {
      cn += " rxChange";
    }

    return String.format(sf_drugNameTemplate, cn, name);
  }
}
