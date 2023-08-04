package org.pharmgkb.pharmcat.reporter.format.html;

import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;

import static org.pharmgkb.pharmcat.reporter.caller.DpydCaller.isDpyd;


/**
 * This class represents a row in the recommendations table.
 *
 * @author Mark Woon
 */
public class Report {
  private final String m_id;
  private final DataSource m_source;
  private final boolean m_matched;
  private final List<String> m_urls;
  private SortedSet<GuidelineReport> m_guidelines;
  private boolean m_notCalled;
  private String m_uncalledGenes;
  private SortedSet<Diplotype> m_unmatchedDiplotypes;
  private boolean m_unmatchedInferred;
  private boolean m_unmatchedDpydInferred;


  public Report(DataSource source, DrugReport drugReport) {
    m_id = drugReport.getId();
    m_source = source;
    m_matched = drugReport.isMatched();
    m_urls = drugReport.getUrls();
    m_guidelines = drugReport.getGuidelines();

    if (!m_matched) {
      m_notCalled = drugReport.getGuidelines().stream().noneMatch(GuidelineReport::isReportable);
      if (m_notCalled) {
        m_uncalledGenes = drugReport.getGuidelines().stream()
            .flatMap(guidelineReport ->  guidelineReport.getGeneReports().stream()
                .filter(geneReport -> !geneReport.isReportable()))
            .map(GeneReport::getGeneDisplay)
            .sorted()
            .collect(Collectors.joining(", "));
      } else {
        m_unmatchedDiplotypes = new TreeSet<>();
        for (GuidelineReport guideline : drugReport.getGuidelines()) {
          for (GeneReport geneReport : guideline.getGeneReports()) {
            for (Diplotype dip : geneReport.getRecommendationDiplotypes()) {
              if (dip.isInferred()) {
                if (isDpyd(geneReport)) {
                  m_unmatchedDpydInferred = true;
                } else {
                  m_unmatchedInferred = true;
                }
              }
              m_unmatchedDiplotypes.add(dip);
            }
          }
        }
      }
    }
  }


  public String getId() {
    return m_id;
  }

  public DataSource getSource() {
    return m_source;
  }

  public List<String> getUrls() {
    return m_urls;
  }


  public boolean isMatched() {
    return m_matched;
  }

  public @Nullable SortedSet<GuidelineReport> getGuidelines() {
    return m_guidelines;
  }


  public boolean isNotCalled() {
    return m_notCalled;
  }

  public @Nullable String getUncalledGenes() {
    return m_uncalledGenes;
  }

  public @Nullable SortedSet<Diplotype> getUnmatchedDiplotypes() {
    return m_unmatchedDiplotypes;
  }

  public boolean isUnmatchedInferred() {
    return m_unmatchedInferred;
  }

  public boolean isUnmatchedDpydInferred() {
    return m_unmatchedDpydInferred;
  }
}
