package org.pharmgkb.pharmcat.util;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.MatchLogic;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


/**
 * Class for handling the logic to match {@link MessageAnnotation} objects to their applicable report objects.
 *
 * @author Ryan Whaley
 */
public class MessageMatcher {

  private Collection<MessageAnnotation> m_messages;
  private ReportContext m_reportContext;

  public MessageMatcher(Collection<MessageAnnotation> messages, ReportContext reportContext) {
    Preconditions.checkNotNull(messages);
    Preconditions.checkNotNull(reportContext);
    m_messages = ImmutableList.copyOf(messages);
    m_reportContext = reportContext;
  }

  @Nonnull
  public List<MessageAnnotation> match(DrugReport guideline) {
    return m_messages.stream()
        .filter(m -> match(m.getMatches(), guideline))
        .collect(Collectors.toList());
  }


  public boolean match(MatchLogic match, DrugReport report) {

    boolean criteriaPass = !match.getDrugs().isEmpty() && !Collections.disjoint(match.getDrugs(), report.getRelatedDrugs());

    if (criteriaPass && match.getDips().size() > 0) {
      GeneReport geneReport = m_reportContext.getGeneReport(match.getGene());
      criteriaPass = geneReport.getMatcherDiplotypes().size() > 0 &&
          geneReport.getMatcherDiplotypes().stream()
              .map(Diplotype::printBare)
              .anyMatch(b -> match.getDips().contains(b));
    }
    
    if (criteriaPass && match.getVariantsMissing().size() > 0) {
      GeneReport geneReport = m_reportContext.getGeneReport(match.getGene());
      criteriaPass = geneReport.getVariantReports().stream()
          .filter(VariantReport::isMissing)
          .map(VariantReport::getDbSnpId)
          .collect(Collectors.toSet())
          .containsAll(match.getVariantsMissing());    
    }

    return criteriaPass;
  }
}
