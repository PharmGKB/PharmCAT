package org.pharmgkb.pharmcat.reporter.format;

import java.io.IOException;
import java.nio.file.Path;
import org.pharmgkb.pharmcat.reporter.ReportContext;


public abstract class AbstractFormat {
  private final Path f_outputPath;
  private final String f_title;

  public AbstractFormat(Path outputPath, String title) {
    f_outputPath = outputPath;
    f_title = title;
  }

  public abstract void write(ReportContext reportContext) throws IOException;

  public Path getOutputPath() {
    return f_outputPath;
  }

  public String getTitle() {
    return f_title;
  }
}
