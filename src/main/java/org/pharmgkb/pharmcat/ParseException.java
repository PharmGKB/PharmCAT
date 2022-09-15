package org.pharmgkb.pharmcat;

/**
 * This exception is thrown when there's a problem parsing.
 *
 * @author Mark Woon
 */
public class ParseException extends RuntimeException {
  private String m_additionalInfo;


  public ParseException(String msg) {
    super(msg);
  }

  public ParseException(String msg, Exception ex) {
    super(msg, ex);
  }


  public void setAdditionalInfo(String additionalInfo) {
    m_additionalInfo = additionalInfo;
  }

  public String getAdditionalInfo() {
    return m_additionalInfo;
  }

  @Override
  public String toString() {
    if (m_additionalInfo != null) {
      return getMessage() + " (" + m_additionalInfo + ")";
    }
    return getMessage();
  }
}
