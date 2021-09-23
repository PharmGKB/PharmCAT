package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;
import com.google.gson.Gson;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.MatchLogic;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


/**
 * Wrapper class that will load message data and then help match the messages to applicable models
 */
public class MessageList {
  private static final String sf_messagesFile   = "org/pharmgkb/pharmcat/definition/messages.json";
  private final List<MessageAnnotation> f_messages;

  /**
   * Public constructor. Will load message data from the file system
   * @throws IOException can occur when reading the messages file
   */
  public MessageList() throws IOException {
    try (BufferedReader reader = Files.newBufferedReader(PathUtils.getPathToResource(sf_messagesFile))) {
      MessageAnnotation[] messages = new Gson().fromJson(reader, MessageAnnotation[].class);
      f_messages = Arrays.asList(messages);
    }
  }

  /**
   * This method will go through all messages and add any matching {@link MessageAnnotation} objects to the
   * {@link GeneReport}
   * @param gene the {@link GeneReport} to possibly add messages to
   */
  public void addMatchingMessagesTo(GeneReport gene) {
    if (!gene.isCalled()) {
      // puposely don't apply any messages if the gene is not called
      return;
    }

    List<MessageAnnotation> messages = f_messages.stream()
        .filter(m -> match(m, gene))
        .collect(Collectors.toList());
    gene.addMessages(messages);
  }

  /**
   * Method that will determine if a message's {@link MatchLogic} applies to a {@link GeneReport}
   * @param message match logic from a {@link MessageAnnotation}
   * @param gene a {@link GeneReport} to try to match to
   * @return true if the given logic matches
   */
  public static boolean match(MessageAnnotation message, GeneReport gene) {
    MatchLogic match = message.getMatches();

    if (!Objects.equals(match.getGene(), gene.getGene())) {
      return false;
    }

    boolean passHapMatchCriteria = match.getHapsCalled().isEmpty()
        || match.getHapsCalled().stream().allMatch(gene::hasHaplotype);
    boolean passHapMissingCriteria = match.getHapsMissing().isEmpty()
        || gene.getUncalledHaplotypes().containsAll(match.getHapsMissing());
    boolean passVariantMatchCriteria = StringUtils.isBlank(match.getVariant())
        || gene.findVariantReport(match.getVariant()).map((v) -> !v.isMissing()).orElse(false);
    boolean passVariantMissingCriteria = match.getVariantsMissing().isEmpty()
        || match.getVariantsMissing().stream().allMatch((r) -> gene.findVariantReport(r).map(VariantReport::isMissing).orElse(false));
    boolean passDipMatchCriteria = match.getDips().isEmpty()
        || match.getDips().stream().allMatch(d -> gene.getMatcherDiplotypes().stream().anyMatch(e -> e.printBare().equals(d)));

    // "ambiguity" messages only apply when the gene unphased, other criteria still apply too
    boolean passAmbiguityCriteria = !message.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY)
        || (!match.getDips().isEmpty() && !gene.isPhased());

    // if it's an "ambiguity" message and a variant is specified then that variant must be het
    boolean passAmbiguityHetCallCriteria = !message.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY)
        || StringUtils.isBlank(match.getVariant())
        || gene.findVariantReport(match.getVariant()).map(VariantReport::isHetCall).orElse(false);

    return passHapMatchCriteria && passHapMissingCriteria && passVariantMatchCriteria && passVariantMissingCriteria
        && passDipMatchCriteria && passAmbiguityCriteria && passAmbiguityHetCallCriteria;
  }

  /**
   * Add all matching messages to the given {@link DrugReport}. Pass in the {@link ReportContext} so we can look up
   * information about related genes.
   * @param drugReport a {@link DrugReport} to add message annotations to
   * @param reportContext the report context to pull related information from
   */
  public void match(DrugReport drugReport, ReportContext reportContext) {
    List<MessageAnnotation> matchedMessages = f_messages.stream()
        .filter(m -> match(m, drugReport, reportContext))
        .collect(Collectors.toList());
    drugReport.addMessages(matchedMessages);
  }

  /**
   * See if the supplied {@link MessageAnnotation} applies to the given {@link DrugReport}.
   *
   * <strong>NOTE:</strong> This method assumes that {@link MessageAnnotation} objects have already been assigned to
   * {@link GeneReport} objects.
   * @param message a message annotation to test for a match
   * @param report a drug report to possibly add the message annotation to
   * @param reportContext the report context to look up related gene information
   * @return true if the message is a match, false otherwise
   */
  private boolean match(MessageAnnotation message, DrugReport report, ReportContext reportContext) {
    MatchLogic match = message.getMatches();

    // if there's no drug specified don't continue
    if (match.getDrugs().isEmpty()) return false;

    // if the message doesn't mention a drug in this message annotation don't continue
    if (Collections.disjoint(match.getDrugs(), report.getRelatedDrugs())) {
      return false;
    }
    // at this point we know the drug is a match, now check the other criteria

    if (StringUtils.isBlank(match.getGene())) return true;

    GeneReport geneReport = reportContext.findGeneReport(match.getGene())
        .orElseThrow(() -> new RuntimeException("Unexpected missing gene report for [" + match.getGene() + "]"));

    return geneReport.getMessages().stream()
        .map(MessageAnnotation::getName)
        .anyMatch((n) -> n.equals(message.getName()));
  }
}
