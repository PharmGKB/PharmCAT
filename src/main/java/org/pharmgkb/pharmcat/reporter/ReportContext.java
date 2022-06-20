package org.pharmgkb.pharmcat.reporter;

import java.io.IOException;
import java.util.Collection;
import java.util.Date;
import java.util.List;
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nonnull;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.definition.MessageList;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.cpic.Drug;
import org.pharmgkb.pharmcat.reporter.model.result.DrugReport;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.reporter.model.result.Genotype;
import org.pharmgkb.pharmcat.util.CliUtils;


/**
 * This class acts as a central context for all data needed to generate the final report.
 *
 * It currently gathers
 * <ul>
 *   <li>{@link GeneCall} objects from the named allele matcher</li>
 *   <li>{@link DrugReport} objects from dosing guideline annotations</li>
 *   <li>Allele definitions on a per-gene basis</li>
 * </ul>
 *
 * @author greytwist
 * @author Ryan Whaley
 */
public class ReportContext {
  @Expose
  @SerializedName("title")
  private final String f_title;
  @Expose
  @SerializedName("generatedOn")
  private final Date f_generatedOn = new Date();
  @Expose
  @SerializedName("pharmcatVersion")
  private final String f_pharmcatVersion = CliUtils.getVersion();
  @Expose
  @SerializedName("genes")
  private final SortedSet<GeneReport> f_geneReports;
  @Expose
  @SerializedName("drugs")
  private final SortedSet<DrugReport> f_drugReports = new TreeSet<>();

  /**
   * Public constructor. Compiles all the incoming data into useful objects to be held for later reporting
   * @param geneReports {@link GeneReport} objects, non-null but can be empty
   */
  public ReportContext(Collection<GeneReport> geneReports, String title) throws IOException {
    f_title = title;

    MessageList messageList = new MessageList();

    // add GeneReports from the Phenotyper
    f_geneReports = new TreeSet<>(geneReports);

    // get CPIC drug data
    DrugCollection drugCollection = new DrugCollection();
    // get DPWG/PharmGKB drug data
    PgkbGuidelineCollection pgkbGuidelineCollection = new PgkbGuidelineCollection();

    // go through all CPIC drugs
    for (Drug drug : drugCollection.listReportable()) {
      DrugReport drugReport = createOrFindDrugReport(drug);
      drugReport.addDrugData(drug, this);

      // add matching recommendations
      List<Genotype> possibleGenotypes = makePossibleGenotypes(drugReport.getRelatedGeneSymbols());
      for (Genotype genotype : possibleGenotypes) {
        drugReport.matchAnnotationsToGenotype(genotype, drug);
      }
      f_drugReports.add(drugReport);

      // add the inverse relationship to gene reports
      for (String gene : drugReport.getRelatedGeneSymbols()) {
        getGeneReport(gene).addRelatedDrugs(drugReport);
      }
    }

    // go through all DPWG-PharmGKB drugs, we iterate this way because one guideline may have multiple chemicals/drugs
    for (String drugName : pgkbGuidelineCollection.getChemicals()) {
      pgkbGuidelineCollection.findGuidelinePackages(drugName).forEach(guidelinePackage -> {
        DrugReport drugReport = createOrFindDrugReport(drugName);

        // add matching groups for possible genotypes
        List<Genotype> possibleGenotypes = makePossibleGenotypes(guidelinePackage.getGenes());
        for (Genotype genotype : possibleGenotypes) {
          guidelinePackage.match(genotype);
        }

        drugReport.addDrugData(guidelinePackage, this);
        f_drugReports.add(drugReport);
        for (String gene : drugReport.getRelatedGeneSymbols()) {
          getGeneReport(gene).addRelatedDrugs(drugReport);
        }
      });
    }

    // now that all reports are generated, apply the applicable messages
    for (DrugReport drugReport : getDrugReports()) {
      messageList.match(drugReport, this);

      // add message to drug when a related gene has a *1 allele
      boolean hasStarOne = drugReport.getRelatedGeneSymbols().stream()
          .flatMap((s) -> getGeneReport(s).getReporterDiplotypes().stream())
          .anyMatch((d) -> d.hasAllele("*1"));
      if (hasStarOne) {
        drugReport.addMessage(new MessageAnnotation(
            MessageAnnotation.TYPE_NOTE,
            "The *1 allele assignment is characterized by the absence of variants that are included in the " +
                "underlying allele definitions by either position being reference or missing."
        ));
      }

      // add a message for any gene that has missing data
      drugReport.getRelatedGeneSymbols().stream()
          .filter((s) -> !getGeneReport(s).isOutsideCall() && getGeneReport(s).isMissingVariants().equals(GeneReport.YES))
          .map((s) -> "Some position data used to define " + s + " alleles is missing which may change the matched " +
              "genotype. See the gene section for " + s + " for more information.")
          .forEach((m) -> drugReport.addMessage(new MessageAnnotation(MessageAnnotation.TYPE_NOTE, m)));
    }
  }

  public SortedSet<DrugReport> getDrugReports() {
    return f_drugReports;
  }

  public Collection<GeneReport> getGeneReports() {
    return f_geneReports;
  }

  private DrugReport createOrFindDrugReport(@Nonnull Drug drug) {
    return createOrFindDrugReport(drug.getDrugName());
  }

  private DrugReport createOrFindDrugReport(@Nonnull String drugName) {
    Preconditions.checkNotNull(drugName);
    return getDrugReports().stream()
        .filter(r -> r.getName().equalsIgnoreCase(drugName))
        .findFirst()
        .orElse(new DrugReport(drugName));
  }

  /**
   * Find a {@link GeneReport} based on the gene symbol
   * @param geneSymbol a gene symbol
   */
  public Optional<GeneReport> findGeneReport(String geneSymbol) {
    return getGeneReports().stream().filter(r -> r.getGene().equals(geneSymbol)).findFirst();
  }

  /**
   * Finds the {@link GeneReport} for the given gene symbol and will throw a RuntimeException if it's not found,
   * effectively guaranteeing a non-null result
   * @param geneSymbol a gene symbol to find a report for
   * @return a GeneReport object
   * @throws RuntimeException if the desired gene report does not exist
   */
  public GeneReport getGeneReport(String geneSymbol) {
    return findGeneReport(geneSymbol)
        .orElseThrow(() -> new RuntimeException("No gene exists for " + geneSymbol));
  }

  /**
   * Find a {@link DrugReport} for the drug with the given name.
   * @param drugName the name of the drug to find a report for
   * @return an Optional {@link DrugReport}
   */
  private Optional<DrugReport> findDrugReport(String drugName) {
    return f_drugReports.stream()
        .filter(r -> r.getRelatedDrugs().contains(drugName))
        .findFirst();
  }

  /**
   * Gets a {@link DrugReport} for the drug with the given name. Will throw a {@link RuntimeException} if the drug is
   * not found.
   * @param drugName the name of the drug to find
   * @return a non-null {@link DrugReport}
   */
  public DrugReport getDrugReport(String drugName) {
    return findDrugReport(drugName).orElseThrow(() -> new RuntimeException("No drug exists for " + drugName));
  }

  private List<Genotype> makePossibleGenotypes(Collection<String> geneSymbols) {
    List<GeneReport> geneReports = getGeneReports().stream()
        .filter(r -> geneSymbols.contains(r.getGene()))
        .toList();
    return Genotype.makeGenotypes(geneReports);
  }

  public String getTitle() {
    return f_title;
  }

  public Date getGeneratedOn() {
    return f_generatedOn;
  }

  public String getPharmcatVersion() {
    return f_pharmcatVersion;
  }
}
