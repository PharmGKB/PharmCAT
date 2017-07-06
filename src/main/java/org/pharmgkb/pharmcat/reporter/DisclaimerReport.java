package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import org.pharmgkb.common.io.util.CliHelper;


/**
 * CLI to write out a single HTML version of the disclaimer text
 *
 * @author Ryan Whaley
 */
public class DisclaimerReport {

  public static void main(String[] args) {
    CliHelper cliHelper = new CliHelper(MethodHandles.lookup().lookupClass())
        .addOption("o", "output-file", "path to write disclaimer file to", true, "o");

    if (!cliHelper.parse(args)) {
      System.exit(1);
    }

    try {
      HtmlReportGenerator.writeDisclaimerReport(cliHelper.getPath("o"));

      System.out.println("wrote to " + cliHelper.getPath("o"));
    }
    catch (IOException e) {
      e.printStackTrace();
      System.exit(1);
    }
  }
}
