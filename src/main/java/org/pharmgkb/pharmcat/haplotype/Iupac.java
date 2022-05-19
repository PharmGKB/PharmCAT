package org.pharmgkb.pharmcat.haplotype;

import java.util.List;
import java.util.Objects;
import java.util.regex.Pattern;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * Enumeration of IUPAC DNA codes.
 *
 * @author Mark Woon
 */
public enum Iupac {
  A("A", "A", false),
  C("C", "C", false),
  G("G", "G", false),
  T("T", "T", false),
  R("R", "[AG]", true),
  Y("Y", "[CT]", true),
  S("S", "[GC]", true),
  W("W", "[AT]", true),
  K("K", "[GT]", true),
  M("M", "[AC]", true),
  B("B", "[CGT]", true),
  D("D", "[AGT]", true),
  H("H", "[ACT]", true),
  V("V", "[ACG]", true),
  N("N", "[ACGT]", true),
  DEL("-", "del", false);


  private final String m_code;
  private Pattern m_pattern;
  private final boolean m_ambiguity;


  Iupac(String code, @Nullable String regex, boolean isAmbiguity) {
    m_code = code;
    m_ambiguity = isAmbiguity;
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

  public boolean isAmbiguity() {
    return m_ambiguity;
  }

  public List<String> getBases() {
    if (this == DEL) {
      return ImmutableList.of();
    } else if (getRegex().length() == 1) {
      return ImmutableList.of(getRegex());
    } else {
      return ImmutableList.copyOf(getRegex().substring(1, getRegex().length() - 1).split("(?!^)"));
    }
  }

  public static Iupac lookup(String value) {
    Preconditions.checkNotNull(value);

    if (value.equals("-")) {
      return DEL;
    }
    return valueOf(value.toUpperCase());
  }

  public static boolean isWobble(String allele) {
    if (allele.length() == 1) {
      return Objects.requireNonNull(lookup(allele)).isAmbiguity();
    }
    return false;
  }
}
