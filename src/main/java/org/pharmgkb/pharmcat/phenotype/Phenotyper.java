package org.pharmgkb.pharmcat.phenotype;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.pharmgkb.pharmcat.definition.MessageList;
import org.pharmgkb.pharmcat.definition.PhenotypeMap;
import org.pharmgkb.pharmcat.definition.ReferenceAlleleMap;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.DrugCollection;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class takes genotyping information from the {@link NamedAlleleMatcher} and from outside allele calls then
 * interprets them into phenotype assignments, diplotype calls, and other information needed for subsequent use in the
 * {@link ReportContext}. The data is compiled into {@link GeneReport} objects which can then serialized
 */
public class Phenotyper {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private final SortedSet<GeneReport> f_geneReports = new TreeSet<>();


  /**
   * Public constructor. This needs {@link GeneCall} objects from the {@link NamedAlleleMatcher} and {@link OutsideCall}
   * objects coming from other allele calling sources. This relies on reading definition files as well.
   * @param geneCalls a List of {@link GeneCall} objects
   * @param outsideCalls a List of {@link OutsideCall} objects
   */
  public Phenotyper(List<GeneCall> geneCalls, List<OutsideCall> outsideCalls,
      @Nullable Map<String, Collection<String>> variantWarnings) {
    ReferenceAlleleMap referenceAlleleMap = new ReferenceAlleleMap();
    PhenotypeMap phenotypeMap = new PhenotypeMap();

    for (GeneCall geneCall : geneCalls) {
      GeneReport geneReport = new GeneReport(geneCall);
      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
          geneReport.getGene(),
          phenotypeMap.lookup(geneReport.getGene()).orElse(null),
          referenceAlleleMap.get(geneReport.getGene()));
      geneReport.setDiplotypes(diplotypeFactory, geneCall);

      f_geneReports.add(geneReport);
    }

    for (OutsideCall outsideCall : outsideCalls) {
      GeneReport geneReport = findGeneReport(outsideCall.getGene())
          .orElse(null);

      if (geneReport != null && geneReport.getCallSource() != CallSource.OUTSIDE) {
        String matcherCall = String.join("; ", geneReport.printDisplayCalls());

        // remove the existing call so we can use the new outside call
        removeGeneReport(outsideCall.getGene());

        // use this new outside call instead
        geneReport = new GeneReport(outsideCall);
        f_geneReports.add(geneReport);

        // warn the user of the conflict
        geneReport.addMessage(new MessageAnnotation(
            MessageAnnotation.TYPE_NOTE,
            "prefer-sample-data",
            "VCF data would call " + matcherCall + " but it has been ignored in favor of an outside call. If you want to use the matcher call for this gene then remove the gene from the outside call data."
        ));
      }

      if (geneReport == null) {
        geneReport = new GeneReport(outsideCall);
        f_geneReports.add(geneReport);
      }

      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(
          geneReport.getGene(),
          phenotypeMap.lookup(geneReport.getGene()).orElse(null),
          referenceAlleleMap.get(geneReport.getGene()));
      geneReport.setDiplotypes(diplotypeFactory, outsideCall);

      f_geneReports.add(geneReport);
    }

    Set<String> unspecifiedGenes = listUnspecifiedGenes();
    for (String geneSymbol : unspecifiedGenes) {
      GeneReport geneReport = new GeneReport(geneSymbol);
      DiplotypeFactory diplotypeFactory = new DiplotypeFactory(geneSymbol, phenotypeMap.lookup(geneSymbol).orElse(null), null);
      geneReport.setUnknownDiplotype(diplotypeFactory);
      f_geneReports.add(geneReport);
    }

    f_geneReports.forEach(r -> r.addVariantWarningMessages(variantWarnings));

    try {
      MessageList messageList = new MessageList();
      getGeneReports().forEach(messageList::addMatchingMessagesTo);
    } catch (IOException ex) {
      throw new RuntimeException("Could not load messages", ex);
    }
  }

  /**
   * Write the collection of {@link GeneReport} objects
   * @param outputPath the path to write a JSON file of data to
   * @throws IOException can occur from writing the file to the filesystem
   */
  public void write(Path outputPath) throws IOException {
    try (BufferedWriter writer = Files.newBufferedWriter(outputPath, StandardCharsets.UTF_8)) {
      Gson gson = new GsonBuilder().serializeNulls().excludeFieldsWithoutExposeAnnotation()
          .setPrettyPrinting().create();
      writer.write(gson.toJson(f_geneReports));
      sf_logger.info("Writing Phenotyper JSON to " + outputPath);
    }
  }

  /**
   * Get the {@link GeneReport} objects that were created from call data in the constructor
   * @return a SortedSet of {@link GeneReport} objects
   */
  public SortedSet<GeneReport> getGeneReports() {
    return f_geneReports;
  }

  /**
   * Find a {@link GeneReport} based on the gene symbol
   * @param geneSymbol a gene symbol
   */
  public Optional<GeneReport> findGeneReport(String geneSymbol) {
    return getGeneReports().stream().filter(r -> r.getGene().equals(geneSymbol)).findFirst();
  }

  private void removeGeneReport(String geneSymbol) {
    findGeneReport(geneSymbol).ifPresent(f_geneReports::remove);
  }

  /**
   * Read gene resport information from a given JSON file path. This should be the JSON output of this class.
   * @param filePath a path to an existing JSON file
   * @return a List of {@link GeneReport} objects
   * @throws IOException can occur if file is unable to be read
   */
  public static List<GeneReport> readGeneReports(Path filePath) throws IOException {
    Preconditions.checkNotNull(filePath);
    Preconditions.checkArgument(filePath.toFile().exists());
    Preconditions.checkArgument(filePath.toFile().isFile());
    Gson gson = new GsonBuilder().create();
    try (BufferedReader reader = Files.newBufferedReader(filePath)) {
      return Arrays.asList(gson.fromJson(reader, GeneReport[].class));
    }
  }

  public Set<String> listUnspecifiedGenes() {
    try {
      DrugCollection cpicDrugs = new DrugCollection();
      PgkbGuidelineCollection dpwgDrugs = new PgkbGuidelineCollection();

      Set<String> unspecifiedGenes = new HashSet<>();
      unspecifiedGenes.addAll(cpicDrugs.getAllReportableGenes());
      unspecifiedGenes.addAll(dpwgDrugs.getGenes());
      f_geneReports.stream()
          .map(GeneReport::getGene)
          .forEach(unspecifiedGenes::remove);
      return unspecifiedGenes;
    } catch (IOException ex) {
      throw new RuntimeException("Error reading drug data", ex);
    }
  }
}
