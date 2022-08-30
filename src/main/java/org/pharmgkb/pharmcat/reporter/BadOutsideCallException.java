package org.pharmgkb.pharmcat.reporter;

/**
 * This exception indicates a prolem with an outside call.
 *
 * @author Mark Woon
 */
public class BadOutsideCallException extends RuntimeException {
  public BadOutsideCallException(String msg) {
    super(msg);
  }
}
