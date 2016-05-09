package org.pharmgkb.pharmcat.haplotype;

import java.util.regex.Pattern;
import javax.annotation.Nonnull;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;

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
  R("R", "[AG]"),
  Y("Y", "[CT]"),
  S("S", "[GC]"),
  W("W", "[AT]"),
  K("K", "[GT]"),
  M("M", "[AC]"),
  B("B", "[CGT]"),
  D("D", "[AGT]"),
  H("H", "[ACT]"),
  V("V", "[ACG]"),
  N("N", "[ACGT]"),
  DEL("-", "del");


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

    if (value.equals("-")) {
      return DEL;
    }
    return valueOf(value.toUpperCase());
  }
}
