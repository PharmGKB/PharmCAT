package org.pharmgkb.pharmcat.reporter;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Objects;
import java.util.Optional;
import java.util.stream.Stream;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MatchLogic;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.AnnotationGroup;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.GuidelineReport;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;
import org.pharmgkb.pharmcat.util.DataSerializer;


/**
 * Helper class to help assign custom messages to {@link GeneReport} and {@link DrugReport}.
 */
public class MessageHelper {
  public static final String MESSAGES_JSON_FILE_NAME = "messages.json";
  private static final String sf_messagesFile   = "org/pharmgkb/pharmcat/reporter/" + MESSAGES_JSON_FILE_NAME;
  private final Multimap<String, MessageAnnotation> m_geneMap = HashMultimap.create();
  private final Multimap<String, MessageAnnotation> m_drugMap = HashMultimap.create();


  /**
   * Public constructor. Will load message data from the file system.
   *
   * @throws IOException can occur when reading the messages file
   */
  public MessageHelper() throws IOException {
    try (BufferedReader reader = Files.newBufferedReader(PathUtils.getPathToResource(sf_messagesFile))) {
      MessageAnnotation[] messages = DataSerializer.GSON.fromJson(reader, MessageAnnotation[].class);
      for (MessageAnnotation msg : messages) {
        if (msg.getMatches().getGene() != null) {
          m_geneMap.put(msg.getMatches().getGene(), msg);
        }
        msg.getMatches().getDrugs()
            .forEach((d) -> m_drugMap.put(d, msg));
      }
    }
  }

  /**
   * This method will go through all messages and add any matching {@link MessageAnnotation} objects to the
   * {@link GeneReport}.
   *
   * @param report the {@link GeneReport} to possibly add messages to
   */
  public void addMatchingMessagesTo(GeneReport report) {
    if (!report.isCalled()) {
      // puposely don't apply any messages if the gene is not called
      return;
    }
    if (report.getCallSource() != CallSource.MATCHER) {
      return;
    }
    m_geneMap.get(report.getGene()).stream()
        .filter(m -> matchesGeneReport(m, report))
        .forEach(report::addMessage);
  }

  /**
   * Method that will determine if a message's {@link MatchLogic} applies to a {@link GeneReport}
   * @param message match logic from a {@link MessageAnnotation}
   * @param gene a {@link GeneReport} to try to match to
   * @return true if the given logic matches
   */
  private static boolean matchesGeneReport(MessageAnnotation message, GeneReport gene) {
    MatchLogic match = message.getMatches();

    if (!Objects.equals(match.getGene(), gene.getGene())) {
      return false;
    }

    boolean passHapMatchCriteria = match.getHapsCalled().isEmpty() ||
        match.getHapsCalled().stream().allMatch(gene::hasHaplotype);
    boolean passHapMissingCriteria = match.getHapsMissing().isEmpty() ||
        gene.getUncalledHaplotypes().containsAll(match.getHapsMissing());
    boolean passVariantMatchCriteria = StringUtils.isBlank(match.getVariant()) ||
        gene.findVariantReport(match.getVariant()).map((v) -> !v.isMissing()).orElse(false);
    boolean passVariantMissingCriteria = match.getVariantsMissing().isEmpty() ||
        match.getVariantsMissing().stream()
            .allMatch((r) -> gene.findVariantReport(r).map(VariantReport::isMissing).orElse(false));
    boolean passDipMatchCriteria = match.getDips().isEmpty() ||
        match.getDips().stream()
            .allMatch(d -> gene.getSourceDiplotypes().stream().anyMatch(e -> e.printBare().equals(d)));

    boolean passAmbiguityCriteria = true;
    if (message.getExceptionType().equals(MessageAnnotation.TYPE_AMBIGUITY)) {
      // ambiguity messages with diplotypes only apply if gene is unphased
      if (!match.getDips().isEmpty() && gene.isPhased()) {
        passAmbiguityCriteria = false;
      }
      // ambiguity messages with a variant only apply when that variant is het
      else if (!StringUtils.isBlank(match.getVariant()) &&
          !gene.findVariantReport(match.getVariant()).map(VariantReport::isHetCall).orElse(false)) {
        passAmbiguityCriteria = false;
      }
    }

    return passHapMatchCriteria && passHapMissingCriteria && passVariantMatchCriteria && passVariantMissingCriteria
        && passDipMatchCriteria && passAmbiguityCriteria;
  }

  /**
   * Add all matching messages to the given {@link DrugReport}. Pass in the {@link ReportContext} so we can look up
   * information about related genes.
   *
   * @param drugReport a {@link DrugReport} to add message annotations to
   * @param reportContext the report context to pull related information from
   */
  public void addMatchingMessagesTo(DrugReport drugReport, ReportContext reportContext, DataSource source) {
    Collection<MessageAnnotation> allMessages = m_drugMap.get(drugReport.getName()).stream()
        .filter((ma) -> allowedForSource(ma, source))
        .toList();
    if (allMessages.size() == 0) {
      return;
    }
    List<MessageAnnotation> reportAsGenotype = new ArrayList<>();
    for (MessageAnnotation messageAnnotation : allMessages) {
      if (messageAnnotation.getExceptionType().equals(MessageAnnotation.TYPE_REPORT_AS_GENOTYPE)) {
        reportAsGenotype.add(messageAnnotation);
      } else {
        if (matchDrugReport(messageAnnotation, reportContext, source)) {
          drugReport.addMessage(messageAnnotation);
        }
      }
    }

    if (reportAsGenotype.size() > 0) {
      for (MessageAnnotation msgAnn : reportAsGenotype) {
        String geneSymbol = msgAnn.getMatches().getGene();
        String genotype = null;
        for (GuidelineReport guidelineReport : drugReport.getGuidelines()) {
          if (geneSymbol == null || guidelineReport.getRelatedGenes().contains(geneSymbol)) {
            for (AnnotationGroup annotationGroup : guidelineReport.getAnnotationGroups()) {
              if (genotype == null) {
                genotype = computeGenotype(msgAnn, reportContext, source);
              }
              annotationGroup.addHighlightedVariant(genotype);
            }
          }
        }
      }
    }
  }

  private boolean allowedForSource(MessageAnnotation messageAnnotation, DataSource source) {
    String key = messageAnnotation.getName();
    if (key.contains("cpic-") && source != DataSource.CPIC) {
      return false;
    }
    if (key.contains("dpwg-") && source != DataSource.DPWG) {
      return false;
    }
    return true;
  }


  private String computeGenotype(MessageAnnotation msgAnn, ReportContext reportContext, DataSource source) {
    String geneSymbol = Objects.requireNonNull(msgAnn.getMatches().getGene());
    String rsid = Objects.requireNonNull(msgAnn.getMatches().getVariant());

    GeneReport gr = reportContext.getGeneReport(source, geneSymbol);
    Optional<String> call = Stream.concat(gr.getVariantReports().stream(), gr.getVariantOfInterestReports().stream())
        .filter(v -> v.getDbSnpId() != null && v.getDbSnpId().matches(rsid) && !v.isMissing())
        .map(VariantReport::getCall)
        .findFirst();
    if (call.isEmpty() || StringUtils.isBlank(call.get())) {
      return rsid + ": " + Haplotype.UNKNOWN;
    }
    else {
      return rsid + ": " + call.get().replaceAll("\\|", "/");
    }
  }


  /**
   * See if the supplied {@link MessageAnnotation} applies to the given {@link DrugReport} based on gene-related
   * matches.
   * <p>
   * <strong>NOTE:</strong> This method assumes that {@link MessageAnnotation} objects have already been assigned to
   * {@link GeneReport} objects.
   *
   * @param message a message annotation to test for a match
   * @param reportContext the report context to look up related gene information
   * @return true if the message is a match, false otherwise
   */
  private boolean matchDrugReport(MessageAnnotation message, ReportContext reportContext,
      DataSource source) {
    String gene = message.getMatches().getGene();
    return StringUtils.isBlank(gene) || reportContext.getGeneReport(source, gene)
        .hasMessage(message.getName());
  }
}
