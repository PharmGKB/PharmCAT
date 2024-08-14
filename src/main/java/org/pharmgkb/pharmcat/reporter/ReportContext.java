package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.PrescribingGuidanceSource;
import org.pharmgkb.pharmcat.reporter.model.pgkb.GuidelinePackage;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.CliUtils;


/**
 * This class acts as a central context for all data needed to generate the final report.
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class ReportContext {
  @Expose
  @SerializedName("title")
  private final String f_title;
  @Expose
  @SerializedName("timestamp")
  private final Date m_timestamp = new Date();
  @Expose
  @SerializedName("pharmcatVersion")
  private final String f_pharmcatVersion = CliUtils.getVersion();
  @Expose
  @SerializedName("dataVersion")
  private String m_dataVersion;
  @Expose
  @SerializedName("genes")
  private final SortedMap<DataSource, SortedMap<String, GeneReport>> m_geneReports;
  @Expose
  @SerializedName("drugs")
  private final SortedMap<PrescribingGuidanceSource, SortedMap<String, DrugReport>> m_drugReports = new TreeMap<>();
  @Expose
  @SerializedName("messages")
  private final List<MessageAnnotation> f_messages = new ArrayList<>();

  /**
   * Public constructor. Compiles all the incoming data into useful objects to be held for later reporting.
   *
   * @param geneReports {@link GeneReport} objects, non-null but can be empty
   * @param title the optional text to show as a user-friendly title or identifier for this report
   */
  public ReportContext(Env env, SortedMap<DataSource, SortedMap<String, GeneReport>> geneReports, String title) throws IOException {
    f_title = title;
    m_geneReports = geneReports;

    m_dataVersion = validateVersions(env.getDrugs());

    for (PrescribingGuidanceSource dataSourceType : PrescribingGuidanceSource.values()) {
      Map<String, DrugReport> drugReports = m_drugReports.computeIfAbsent(dataSourceType, (s) -> new TreeMap<>());
      // go through all drugs, we iterate this way because one guideline may have multiple chemicals/drugs
      for (String drugName : env.getDrugs().getGuidelineMap().keys()) {
        List<GuidelinePackage> guidelinePackages = env.getDrugs().findGuidelinePackages(drugName, dataSourceType);
        if (guidelinePackages != null && !guidelinePackages.isEmpty()) {
          DrugReport newDrugReport = new DrugReport(drugName, guidelinePackages, this);
          drugReports.put(drugName.toLowerCase(), newDrugReport);
        }
      }
    }

    // now that all reports are generated, apply applicable messages
    MessageHelper messageHelper = env.getMessageHelper();
    // to gene reports
    geneReports.values().stream()
        .flatMap((m) -> m.values().stream())
        .forEach(messageHelper::addMatchingMessagesTo);
    // to drug reports
    for (PrescribingGuidanceSource source : m_drugReports.keySet()) {
      for (DrugReport drugReport : m_drugReports.get(source).values()) {
        messageHelper.addMatchingMessagesTo(drugReport, this, source);

        // add a message for any gene that has missing data
        drugReport.getRelatedGeneSymbols().stream()
            .map((s) -> getGeneReport(source, s))
            .filter((gr) -> gr != null && !gr.isOutsideCall() && gr.isMissingVariants() && !gr.isNoData())
            .forEach((gr) -> drugReport.addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE,
                "missing-variants",
                "Some position data used to define " + gr.getGeneDisplay() +
                    " alleles is missing which may change the matched genotype. See <a href=\"#" +
                    gr.getGeneDisplay() + "\">" + gr.getGeneDisplay() +
                    "</a> in Section III for for more information.")));
      }
    }
  }

  private String validateVersions(PgkbGuidelineCollection guidelineCollection) {
    Set<String> observedVersions = new HashSet<>();
    List<GeneReport> ungroupedGeneReports = m_geneReports.values().stream().flatMap((m) -> m.values().stream()).toList();

    for (GeneReport geneReport : ungroupedGeneReports) {
        if (geneReport.getAlleleDefinitionVersion() != null) {
          observedVersions.add(geneReport.getAlleleDefinitionVersion());
        }
        if (geneReport.getPhenotypeVersion() != null) {
          observedVersions.add(geneReport.getPhenotypeVersion());
        }
    }

    // check drug data
    if (guidelineCollection.getVersion() != null) {
      observedVersions.add(guidelineCollection.getVersion());
    }

    if (observedVersions.isEmpty()) {
      return TextConstants.NA;
    } else {
      if (observedVersions.size() > 1) {
        addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "multiple-versions",
            "Multiple versions used to generate gene and drug reports: " + observedVersions + "."
        ));
      }
      return String.join(", ", observedVersions);
    }
  }



  /**
   * Gets the set of all {@link DrugReport} objects that hold drug information and their recommendations.
   *
   * @return a map of {@link DrugReport} objects
   */
  public Map<PrescribingGuidanceSource, SortedMap<String, DrugReport>> getDrugReports() {
    return m_drugReports;
  }

  /**
   * Gets the set of all {@link GeneReport} objects that are reported in this context
   */
  public SortedMap<DataSource, SortedMap<String, GeneReport>> getGeneReports() {
    return m_geneReports;
  }

  public List<DrugReport> getDrugReports(String drug) {
    return m_drugReports.keySet().stream()
        .map((k) -> m_drugReports.get(k).get(drug))
        .filter(Objects::nonNull)
        .toList();
  }

  public @Nullable DrugReport getDrugReport(PrescribingGuidanceSource type, String drug) {
    return m_drugReports.get(type).get(drug);
  }

  public List<GeneReport> getGeneReports(String gene) {
    return m_geneReports.keySet().stream()
        .map((k) -> m_geneReports.get(k).get(gene))
        .filter(Objects::nonNull)
        .toList();
  }

  public @Nullable GeneReport getGeneReport(DataSource source, String gene) {
    return m_geneReports.get(source).get(gene);
  }

  public @Nullable GeneReport getGeneReport(PrescribingGuidanceSource source, String gene) {
    return m_geneReports.get(source.getPhenoSource()).get(gene);
  }

  /**
   * The user-friendly title for the report.
   *
   * @return the title string
   */
  public String getTitle() {
    return f_title;
  }

  /**
   * Gets the timestamp this context was compiled.
   *
   * @return the timestamp this context was compiled
   */
  public Date timestamp() {
    return m_timestamp;
  }

  /**
   * Gets the PharmCAT version tag this context was created with
   * @return a version tag string in the form vX.Y
   */
  public String getPharmcatVersion() {
    return f_pharmcatVersion;
  }

  public String getDataVersion() {
    return m_dataVersion;
  }

  public List<MessageAnnotation> getMessages() {
    return f_messages;
  }

  public void addMessage(MessageAnnotation message) {
    f_messages.add(message);
  }
}
