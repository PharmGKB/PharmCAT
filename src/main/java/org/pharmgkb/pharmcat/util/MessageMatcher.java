package org.pharmgkb.pharmcat.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import com.google.common.collect.ImmutableList;
import com.google.gson.Gson;
import org.apache.commons.lang3.StringUtils;
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
        .filter(m -> match(m, guideline, reportContext))
        .collect(Collectors.toList());
    guideline.addMessages(matchedMessages);
  }


  public boolean match(MessageAnnotation message, DrugReport report, ReportContext reportContext) {
    MatchLogic match = message.getMatches();

    boolean criteriaPass = !match.getDrugs().isEmpty() && !Collections.disjoint(match.getDrugs(), report.getRelatedDrugs());

    if (StringUtils.isBlank(match.getGene())) return criteriaPass;

    GeneReport geneReport = reportContext.findGeneReport(match.getGene())
        .orElseThrow(() -> new RuntimeException("Unexpected missing gene report for [" + match.getGene() + "]"));

    // see if there is a matching diplotype call for message's diplotype
    if (criteriaPass && match.getDips().size() > 0) {
      criteriaPass = geneReport.getMatcherDiplotypes() != null && geneReport.getMatcherDiplotypes().size() > 0 &&
          geneReport.getMatcherDiplotypes().stream()
              .map(Diplotype::printBare)
              .anyMatch(b -> match.getDips().contains(b));
    }

    // see if the variant specified in the message is missing
    if (criteriaPass && match.getVariantsMissing().size() > 0) {
      criteriaPass = geneReport.getVariantReports().stream()
          .filter(VariantReport::isMissing)
          .map(VariantReport::getDbSnpId)
          .collect(Collectors.toSet())
          .containsAll(match.getVariantsMissing());    
    }

    // see if there is a heterozygous call for the given RSID in the ambiguity message annotation
    if (criteriaPass && message.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY) && StringUtils.isNotBlank(match.getVariant())) {
      criteriaPass = geneReport.findVariantReport(match.getVariant())
          .map(VariantReport::isHetCall)
          .orElse(false);
    }

    return criteriaPass;
  }
}
