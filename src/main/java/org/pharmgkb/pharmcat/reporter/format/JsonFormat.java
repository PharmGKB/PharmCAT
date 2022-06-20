package org.pharmgkb.pharmcat.reporter.format;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.pharmgkb.pharmcat.reporter.ReportContext;



/**
 * A JSON-formatted version of {@link ReportContext} data.
 */
public class JsonFormat extends AbstractFormat {
  private static final Gson sf_gson = new GsonBuilder().serializeNulls().excludeFieldsWithoutExposeAnnotation()
      .setPrettyPrinting().create();

  public JsonFormat(Path outputPath) {
    super(outputPath);
  }

  @Override
  public void write(ReportContext reportContext) throws IOException {
    try (BufferedWriter writer = Files.newBufferedWriter(getOutputPath(), StandardCharsets.UTF_8)) {
      writer.write(sf_gson.toJson(reportContext));
    }
  }
}
