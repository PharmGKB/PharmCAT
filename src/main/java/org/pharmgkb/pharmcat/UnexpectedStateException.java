package org.pharmgkb.pharmcat;

/**
 * Exception to throw when consuming data that isn't structured as expected
 *
 * @author Ryan Whaley
 */
public class UnexpectedStateException extends RuntimeException {

  public UnexpectedStateException(String msg) {
    super(msg);
  }

}
