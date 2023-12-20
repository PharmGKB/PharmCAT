package org.pharmgkb.pharmcat.reporter.model;

import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.function.Predicate;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;
import org.checkerframework.checker.nullness.qual.NonNull;
import org.pharmgkb.common.util.ComparisonChain;


/**
 * This class is meant to house a message that needs to show on the final PharmCAT report based on some sort of logic
 * to match it to a particular section.
 *
 * @author Lester Carter
 * @author Ryan Whaley
 */
public class MessageAnnotation implements Comparable<MessageAnnotation> {
  public static final String TYPE_AMBIGUITY = "ambiguity";
  public static final String TYPE_COMBO = "combo-partial";
  private static final String TYPE_EXTRA_POSITION = "extra-position-notes";
  private static final String TYPE_FOOTNOTE = "footnote";
  public static final String TYPE_NOTE = "note";
  public static final String TYPE_REPORT_AS_GENOTYPE = "report-as-genotype";
  public static final String TYPE_NONMATCH = "non-match";
  public static Predicate<MessageAnnotation> isFootnote = m -> m.getExceptionType().equals(TYPE_FOOTNOTE);
  public static Predicate<MessageAnnotation> isExtraPositionNote = m -> m.getExceptionType().equals(TYPE_EXTRA_POSITION);
  public static Predicate<MessageAnnotation> isMessage = m -> !m.getExceptionType().equals(TYPE_FOOTNOTE) &&
      !m.getExceptionType().equals(TYPE_EXTRA_POSITION) &&
      !m.getExceptionType().equals(TYPE_REPORT_AS_GENOTYPE) ;

  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults();
  private static final int sf_rowLength = 12;


  /**
   * This constructor parses a TSV row from the
   * <a href="https://docs.google.com/spreadsheets/d/1MkWV6TlJTnw-KRNWeylyUJAUocCgupcJLlmV2fRdtcM">PharmCAT Message
   * Annotations</a> sheet.
   *
   * @param row a TSV row from the messages sheet
   * @throws IllegalArgumentException can occur if the line is not in the expected form
   */
  public MessageAnnotation(String row) throws IllegalArgumentException {
    String[] fields = row.split("\\t");

    if (fields.length < sf_rowLength) {
      throw new IllegalArgumentException("Row not of expected length "+sf_rowLength);
    }

    m_name = fields[0];

    m_matches = new MatchLogic();
    if (fields[1].contains(",")) {
      throw new IllegalArgumentException("Cannot have multiple genes: " + fields[1]);
    }
    m_matches.setGene(fields[1]);
    m_matches.setHapsCalled(parseList(fields[2]));
    m_matches.setHapsMissing(parseList(fields[3]));
    m_matches.setVariant(fields[4]);
    if (fields[4].contains(",")) {
      throw new IllegalArgumentException("Cannot have multiple variants: " + fields[4]);
    }
    m_matches.setVariantsMissing(parseList(fields[5]));
    m_matches.setDips(parseList(fields[6]));
    m_matches.setDrugs(parseList(fields[7]));

    m_exceptionType = fields[9];
    m_message = fields[10];
    m_version = fields[11];
  }

  /**
   * Constructor for dynamic message annotations.
   * <p>
   * This should only be used when the message needs to be generated dynamically and cannot come from the PharmCAT
   * Message Annotations sheet.
   */
  public MessageAnnotation(String type, String name, String message) {
    m_exceptionType = type;
    m_name = name;
    m_message = message;
  }


  @Expose
  @SerializedName("rule_name")
  private String m_name;
  @Expose
  @SerializedName("version")
  private String m_version;
  @Expose
  @SerializedName("matches")
  private MatchLogic m_matches;
  @Expose
  @SerializedName("exception_type")
  private String m_exceptionType;
  @Expose
  @SerializedName("message")
  private String m_message;


  public String getName() {
    return m_name;
  }

  public String getVersion() {
    return m_version;
  }

  public MatchLogic getMatches() {
    return m_matches;
  }

  public String getExceptionType() {
    return m_exceptionType;
  }

  public String getMessage() {
    return m_message;
  }


  @Override
  public String toString() {
    return m_message;
  }

  @Override
  public boolean equals(Object obj) {
    if (this == obj) {
      return true;
    }
    if (!(obj instanceof MessageAnnotation other)) {
      return false;
    }

    return Objects.equals(m_name, other.getName()) &&
        Objects.equals(m_version, other.getVersion()) &&
        Objects.equals(m_exceptionType, other.getExceptionType()) &&
        Objects.equals(m_message, other.getMessage());
  }

  @Override
  public int hashCode() {
    return Objects.hash(m_name, m_version, m_exceptionType, m_matches);
  }


  private static List<String> parseList(String value) {
    if (StringUtils.isBlank(value)) {
      return Collections.emptyList();
    }
    return ImmutableList.copyOf(sf_commaSplitter.split(value));
  }

  @Override
  public int compareTo(@NonNull MessageAnnotation o) {
    if (o == this) {
      return 0;
    }
    return new ComparisonChain()
        .compare(m_name, o.getName())
        .compare(m_version, o.getVersion())
        .compare(m_exceptionType, o.getExceptionType())
        .compare(m_message, o.getMessage())
        .compare(m_matches, o.getMatches())
        .result();
  }
}
