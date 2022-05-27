package org.pharmgkb.pharmcat;

import java.lang.invoke.MethodHandles;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This exception represents an error that should be reported to the user.
 * The stack trace is not important!
 */
public class ReportableException extends Exception {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());

  public ReportableException(String msg) {
    super(msg);
  }

  public ReportableException(String... msgs) {
    super(String.join(System.lineSeparator(), msgs));
  }


  public void logException() {
    sf_logger.warn("ReportableException", this);
  }
}
