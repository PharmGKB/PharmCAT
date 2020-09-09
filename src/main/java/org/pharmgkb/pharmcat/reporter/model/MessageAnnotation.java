package org.pharmgkb.pharmcat.reporter.model;

import java.util.ArrayList;
import java.util.List;
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
  public static final String TYPE_GENOTYPE = "report-genotype";
  private static final int sf_rowLength = 12;

  public static Predicate<MessageAnnotation> isFootnote = m -> m.getExceptionType().equals("footnote");
  public static Predicate<MessageAnnotation> isExtraPositionNote = m -> m.getExceptionType().equals("extra-position-notes");
  public static Predicate<MessageAnnotation> isMessage = m -> !m.getExceptionType().equals("footnote") && !m.getExceptionType().equals("extra-position-notes");

  /**
   * constructor based off of a row of text from a TSV
   * @param row a row of text from the messages TSV file
   * @throws RuntimeException can occur if the line is not in the expected form
   */
  public MessageAnnotation(String row) throws RuntimeException {
    String[] fields = row.split("\\t");

    if (fields.length < sf_rowLength) {
      throw new RuntimeException("Row not of expected length "+sf_rowLength);
    }

    m_name = fields[0];

    m_matches = new MatchLogic();
    m_matches.setGene(fields[1]);
    m_matches.setHapsCalled(parseList(fields[2]));
    m_matches.setHapsMissing(parseList(fields[3]));
    m_matches.setVariant(fields[4]);
    m_matches.setVariantsMissing(parseList(fields[5]));
    m_matches.setDips(parseList(fields[6]));
    m_matches.setDrugs(parseList(fields[7]));

    m_exceptionType = fields[9];
    m_message = fields[10];
    m_version = fields[11];
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

  public void setName(String rule_name) {
    m_name = rule_name;
  }

  public String getVersion() {
    return m_version;
  }

  public void setVersion(String version) {
    m_version = version;
  }

  public MatchLogic getMatches() {
    return m_matches;
  }

  public void setMatches(MatchLogic matches) {
    m_matches = matches;
  }

  public String getExceptionType() {
    return m_exceptionType;
  }

  public void setExceptionType(String exceptionType) {
    m_exceptionType = exceptionType;
  }

  public String getMessage() {
    return m_message;
  }

  public void setMessage(String message) {
    m_message = message;
  }

  @Override
  public String toString() {
    return getName();
  }

  private static List<String> parseList(String value) {
    if (StringUtils.isBlank(value)) {
      return new ArrayList<>();
    }

    return ImmutableList.copyOf(Splitter.on(",").trimResults().split(value));
  }
}
