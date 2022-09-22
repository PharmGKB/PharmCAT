package org.pharmgkb.pharmcat.reporter.model;

import java.util.Collections;
import java.util.List;
import java.util.MissingResourceException;
import java.util.Objects;
import java.util.ResourceBundle;
import java.util.function.Predicate;
import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.apache.commons.lang3.StringUtils;


/**
 * This class is meant to house a message that needs to show on the final PharmCAT report based on some sort of logic
 * to match it to a particular section.
 *
 * @author Lester Carter
 * @author Ryan Whaley
 */
public class MessageAnnotation {
  public static final String TYPE_AMBIGUITY = "ambiguity";
  public static final String TYPE_COMBO = "combo-partial";
  private static final String TYPE_EXTRA_POSITION = "extra-position-notes";
  private static final String TYPE_FOOTNOTE = "footnote";
  public static final String TYPE_NOTE = "note";
  public static final String TYPE_GENOTYPE = "report-genotype";
  public static Predicate<MessageAnnotation> isFootnote = m -> m.getExceptionType().equals(TYPE_FOOTNOTE);
  public static Predicate<MessageAnnotation> isExtraPositionNote = m -> m.getExceptionType().equals(TYPE_EXTRA_POSITION);
  public static Predicate<MessageAnnotation> isMessage = m -> !m.getExceptionType().equals(TYPE_FOOTNOTE) && !m.getExceptionType().equals(TYPE_EXTRA_POSITION);

  public static final String MSG_COMBO_NAMING = "combo-naming";
  public static final String MSG_COMBO_UNPHASED = "combo-unphased";
  public static final String MSG_CYP2D6_GENERAL = "cyp2d6-general";
  public static final String MSG_CYP2D6_MODE = "cyp2d6-mode";
  public static final String MSG_STAR1_ALLELE = "star1-allele";

  private static final ResourceBundle sf_resources = ResourceBundle.getBundle("messages");
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults();
  private static final int sf_rowLength = 12;


  /**
   * Loads static message from messages.properties.
   */
  public static MessageAnnotation loadMessage(String key) {
    String name = key;
    try {
      name = sf_resources.getString(key + "_name");
    } catch (MissingResourceException ex) {
      // ignore
    }
    return new MessageAnnotation(TYPE_NOTE, name, sf_resources.getString(key + "_message"));
  }

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
   * Normal constructor for message annotations.
   * It is left to the caller to determine whether and where the message is meant to apply.
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
    if (m_name != null) {
      return m_name;
    }
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
}
