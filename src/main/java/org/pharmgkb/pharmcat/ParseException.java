package org.pharmgkb.pharmcat;

/**
 * This exception is thrown when there's a problem parsing.
 *
 * @author Mark Woon
 */
public class ParseException extends RuntimeException {


  public ParseException(String msg) {
    super(msg);
  }


  public ParseException(String msg, Exception ex) {
    super(msg, ex);
  }
}
