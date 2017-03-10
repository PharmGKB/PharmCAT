package org.pharmgkb.pharmcat.reporter.model;

import java.util.List;
import com.google.common.collect.ImmutableList;


/**
 * Model object to hold all properties needed to display a link to a drug
 *
 * @author Ryan Whaley
 */
public class DrugLink {

  private static final List<String> sf_highlightDrugs = ImmutableList.of(
      "warfarin",
      "ribavirin",
      "peginterferon alfa-2a",
      "peginterferon alfa-2b"
  );

  private String m_name;
  private String m_guidelineId;
  private String m_rxClass;

  public DrugLink(String name, String guidelineId, boolean rxChange, boolean rxPossibly) {
    m_name = name;
    m_guidelineId = guidelineId;
    if (rxChange) {
      m_rxClass = "rxChange";
    } else if (rxPossibly) {
      m_rxClass = "rxPossibly";
    } else {
      m_rxClass = "";
    }
  }

  public String getName() {
    return m_name;
  }

  public String getGuidelineId() {
    return m_guidelineId;
  }

  public String getRxClass() {
    return m_rxClass;
  }

  public boolean isHighlighted() {
    return sf_highlightDrugs.contains(getName());
  }
}
