package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import com.google.common.collect.ImmutableList;
import com.google.gson.Gson;
import org.pharmgkb.common.util.PathUtils;
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
  private static final String sf_messagesFile   = "org/pharmgkb/pharmcat/definition/messages.json";

  private final List<MessageAnnotation> m_messages;

  public MessageMatcher() throws IOException {
    try (BufferedReader reader = Files.newBufferedReader(PathUtils.getPathToResource(sf_messagesFile))) {
      m_messages = ImmutableList.copyOf(new Gson().fromJson(reader, MessageAnnotation[].class));
    }
  }

  public void match(DrugReport guideline, ReportContext reportContext) {
    List<MessageAnnotation> matchedMessages = m_messages.stream()
        .filter(m -> match(m.getMatches(), guideline, reportContext))
        .collect(Collectors.toList());
    guideline.addMessages(matchedMessages);
  }


  public boolean match(MatchLogic match, DrugReport report, ReportContext reportContext) {

    boolean criteriaPass = !match.getDrugs().isEmpty() && !Collections.disjoint(match.getDrugs(), report.getRelatedDrugs());

    if (criteriaPass && match.getDips().size() > 0) {
      GeneReport geneReport = reportContext.getGeneReport(match.getGene());
      criteriaPass = geneReport.getMatcherDiplotypes() != null && geneReport.getMatcherDiplotypes().size() > 0 &&
          geneReport.getMatcherDiplotypes().stream()
              .map(Diplotype::printBare)
              .anyMatch(b -> match.getDips().contains(b));
    }
    
    if (criteriaPass && match.getVariantsMissing().size() > 0) {
      GeneReport geneReport = reportContext.getGeneReport(match.getGene());
      criteriaPass = geneReport.getVariantReports().stream()
          .filter(VariantReport::isMissing)
          .map(VariantReport::getDbSnpId)
          .collect(Collectors.toSet())
          .containsAll(match.getVariantsMissing());    
    }

    return criteriaPass;
  }
}
