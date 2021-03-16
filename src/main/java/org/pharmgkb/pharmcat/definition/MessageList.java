package org.pharmgkb.pharmcat.definition;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import com.google.gson.Gson;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.model.MatchLogic;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
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
        .filter(m -> match(m.getMatches(), gene))
        .collect(Collectors.toList());
    gene.addMessages(messages);
  }

  /**
   * Method that will determine if a message's {@link MatchLogic} applies to a {@link GeneReport}
   * @param match match logic from a {@link MessageAnnotation}
   * @param gene a {@link GeneReport} to try to match to
   * @return true if the given logic matches
   */
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

}
