package org.cpic.haplotype;

import com.google.common.base.Preconditions;

import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import java.util.regex.Pattern;

/**
 * Enumeration of IUPAC DNA codes.
 *
 * @author Mark Woon
 */
public enum Iupac {
  A("A", "A"),
  C("C", "C"),
  G("G", "G"),
  T("T", "T"),
  R("R", "A|G"),
  Y("Y", "C|T"),
  S("S", "G|C"),
  W("W", "A|T"),
  K("K", "G|T"),
  M("M", "A|C"),
  B("B", "C|G|T"),
  D("D", "A|G|T"),
  H("H", "A|C|T"),
  V("V", "A|C|G"),
  N("N", "A|C|G|T"),
  DEL(".", null);


  private String m_code;
  private Pattern m_pattern;


  Iupac(@Nonnull String code, @Nullable String regex) {
    m_code = code;
    if (regex != null) {
      m_pattern = Pattern.compile(regex);
    }
  }

  public String getCode() {
    return m_code;
  }

  public String getRegex() {
    return m_pattern.pattern();
  }


  public static @Nonnull Iupac lookup(@Nonnull String value) {
    Preconditions.checkNotNull(value);

    if (value.equals(".")) {
      return DEL;
    }
    return valueOf(value.toUpperCase());
  }
}
