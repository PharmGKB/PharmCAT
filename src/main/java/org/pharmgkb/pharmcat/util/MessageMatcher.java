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
  public List<MessageAnnotation> match(GeneReport gene) {
    if (!gene.isCalled()) {
      // puposely don't apply any messages if the gene is not called
      return ImmutableList.of();
    }

    return m_messages.stream()
        .filter(m -> match(m.getMatches(), gene))
        .collect(Collectors.toList());
  }

  @Nonnull
  public List<MessageAnnotation> match(DrugReport guideline) {
    return m_messages.stream()
        .filter(m -> match(m.getMatches(), guideline))
        .collect(Collectors.toList());
  }



  public static boolean match(MatchLogic match, GeneReport gene) {

    boolean criteriaPass = !match.getGene().isEmpty() && match.getGene().equals(gene.getGene());

    if (criteriaPass && !match.getHapsCalled().isEmpty()) {
      criteriaPass = match.getHapsCalled().stream().anyMatch(h -> gene.getMatcherDiplotypes().stream().anyMatch(d -> d.hasAllele(h)));
    }

    if (criteriaPass && !match.getHapsMissing().isEmpty()) {
      criteriaPass = match.getHapsMissing().isEmpty()
          || match.getHapsMissing().stream().allMatch(h -> gene.getUncalledHaplotypes().contains(h));
    }

    if (criteriaPass && !match.getDips().isEmpty()) {
      criteriaPass = match.getDips().stream().allMatch(d -> gene.getMatcherDiplotypes().stream().anyMatch(e -> e.printBare().equals(d)));
    }

    if (criteriaPass && !match.getVariantsMissing().isEmpty()) {
      criteriaPass = match.getVariantsMissing().stream().allMatch(v -> gene.getVariantReports().isEmpty() || gene.getVariantReports().stream()
          .anyMatch(a -> a.getDbSnpId() != null && a.getDbSnpId().equals(v) && a.isMissing()));
    }

    return criteriaPass;
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
