package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Map;
import javax.annotation.Nonnull;
import com.github.jknack.handlebars.Handlebars;
import com.github.jknack.handlebars.Template;
import com.github.jknack.handlebars.helper.StringHelpers;
import com.github.jknack.handlebars.io.ClassPathTemplateLoader;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;


/**
 * Generates HTML reports using Handlebars templating system
 *
 * @author Ryan Whaley
 */
public class HtmlReportGenerator {
  private static final String FINAL_REPORT      = "report";
  private static final String sf_templatePrefix = "/org/pharmgkb/pharmcat/reporter";

  /**
   * Generate a final report for a Map of data
   * @param data a Map of data from the reporter system
   * @param filePath the path to write the report to
   */
  public static void writeFinalReport(@Nonnull Map<String,Object> data, @Nonnull Path filePath) throws IOException {
    Handlebars handlebars = new Handlebars(new ClassPathTemplateLoader(sf_templatePrefix));
    StringHelpers.register(handlebars);
    handlebars.registerHelpers(ReportHelpers.class);
    Template template = handlebars.compile(FINAL_REPORT);

    String html = template.apply(data);
    try (BufferedWriter writer = Files.newBufferedWriter(filePath, StandardCharsets.UTF_8)) {
      writer.write(html);
    }
  }
}
