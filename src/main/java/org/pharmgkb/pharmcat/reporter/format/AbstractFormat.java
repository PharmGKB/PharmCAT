package org.pharmgkb.pharmcat.reporter.format;

import java.io.IOException;
import java.nio.file.Path;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.ReportContext;


/**
 * Class to extend for making report formats
 */
public abstract class AbstractFormat {
  private final Path m_outputPath;
  private final Env m_env;

  /**
   * Constructor. Needs the path to write the output to
   * @param outputPath the path to the file to write to
   */
  public AbstractFormat(Path outputPath, Env env) {
    m_outputPath = outputPath;
    m_env = env;
  }

  /**
   * Write the {@link ReportContext} data out to the path from the constructor.
   * @param reportContext a {@link ReportContext} object with data
   * @throws IOException can occur from disk IO
   */
  public abstract void write(ReportContext reportContext) throws IOException;

  public Path getOutputPath() {
    return m_outputPath;
  }

  public Env getEnv() {
    return m_env;
  }
}
