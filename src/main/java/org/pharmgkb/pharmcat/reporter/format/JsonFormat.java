package org.pharmgkb.pharmcat.reporter.format;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.util.DataSerializer;


/**
 * A JSON-formatted version of {@link ReportContext} data.
 */
public class JsonFormat extends AbstractFormat {

  public JsonFormat(Path outputPath) {
    super(outputPath);
  }

  @Override
  public void write(ReportContext reportContext) throws IOException {
    try (BufferedWriter writer = Files.newBufferedWriter(getOutputPath(), StandardCharsets.UTF_8)) {
      writer.write(DataSerializer.GSON.toJson(reportContext));
    }
  }
}
