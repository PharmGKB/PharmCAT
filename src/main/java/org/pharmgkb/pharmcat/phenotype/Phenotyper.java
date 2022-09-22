package org.pharmgkb.pharmcat.phenotype;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import javax.annotation.Nullable;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.reporter.DrugCollection;
import org.pharmgkb.pharmcat.reporter.MessageHelper;
import org.pharmgkb.pharmcat.reporter.PgkbGuidelineCollection;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
import org.pharmgkb.pharmcat.reporter.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.model.result.CallSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.DataSerializer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * This class takes genotyping information from the {@link NamedAlleleMatcher} and from outside allele calls then
 * interprets them into phenotype assignments, diplotype calls, and other information needed for subsequent use in the
 * {@link ReportContext}. The data is compiled into {@link GeneReport} objects which can then serialized
 */
public class Phenotyper {
  private static final Logger sf_logger = LoggerFactory.getLogger(MethodHandles.lookup().lookupClass());
  private static final List<String> sf_cpicGenesWithNoPhenotype = Lists.newArrayList("CYP4F2", "IFNL3");

  @Expose
  @SerializedName("geneReports")
  private final Map<DataSource, SortedMap<String, GeneReport>> m_geneReports = new HashMap<>();


  /**
   * Public constructor. This needs {@link GeneCall} objects from the {@link NamedAlleleMatcher} and {@link OutsideCall}
   * objects coming from other allele calling sources. This relies on reading definition files as well.
   * @param geneCalls a List of {@link GeneCall} objects
   * @param outsideCalls a List of {@link OutsideCall} objects
   */
  public Phenotyper(Env env, List<GeneCall> geneCalls, List<OutsideCall> outsideCalls,
      @Nullable Map<String, Collection<String>> variantWarnings) {
    initialize(geneCalls, outsideCalls, env, DataSource.CPIC, variantWarnings);
    initialize(geneCalls, outsideCalls, env, DataSource.DPWG, variantWarnings);
  }


  private void initialize(List<GeneCall> geneCalls, List<OutsideCall> outsideCalls, Env env, DataSource source,
      @Nullable Map<String, Collection<String>> variantWarnings) {
    SortedMap<String, GeneReport> reportMap = m_geneReports.computeIfAbsent(source, (s) -> new TreeMap<>());

    // matcher calls
    for (GeneCall geneCall : geneCalls) {
      if (env.getPhenotype(geneCall.getGene(), source) == null) {
        if (!(source == DataSource.CPIC && sf_cpicGenesWithNoPhenotype.contains(geneCall.getGene()))) {
          continue;
        }
      }
      GeneReport geneReport = new GeneReport(geneCall, env, source);
      reportMap.put(geneReport.getGene(), geneReport);
    }

    //  outside calls
    for (OutsideCall outsideCall : outsideCalls) {
      GeneReport geneReport = reportMap.get(outsideCall.getGene());
      MessageAnnotation msgAnnotation = null;
      if (geneReport != null) {
        if (geneReport.getCallSource() != CallSource.OUTSIDE) {
          // outside call trumps matcher
          // warn the user of the conflict
          String matcherCall = geneReport.getSourceDiplotypes().stream()
              .sorted()
              .map(Diplotype::printBare)
              .collect(Collectors.joining("; "));
          msgAnnotation = new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "prefer-sample-data",
              "PharmCAT would have called " + matcherCall + " based on the VCF but it has been ignored in " +
                  "favor of the outside call. If you want to use the PharmCAT call for this gene then remove the " +
                  "gene from the outside call data."
          );

        } else {
          // add alternate outside call
          System.out.println("WARNING: Multiple outside calls for " + outsideCall.getGene());
          geneReport.addOutsideCall(outsideCall, env);
          continue;
        }
      }

      geneReport = new GeneReport(outsideCall, env, source);
      if (msgAnnotation != null) {
        geneReport.addMessage(msgAnnotation);
      }
      reportMap.put(geneReport.getGene(), geneReport);
    }

    // all other genes
    for (String geneSymbol : listUnspecifiedGenes(source)) {
      reportMap.put(geneSymbol, GeneReport.unspecifiedGeneReport(geneSymbol, env, source));
    }

    try {
      MessageHelper messageHelper = new MessageHelper();
      reportMap.values().forEach(r -> {
        r.addVariantWarningMessages(variantWarnings);
        messageHelper.addMatchingMessagesTo(r);
      });
    } catch (IOException ex) {
      throw new RuntimeException("Could not load messages", ex);
    }
  }


  public Map<DataSource, SortedMap<String, GeneReport>> getGeneReports() {
    return m_geneReports;
  }

  /**
   * Find a {@link GeneReport} based on the gene symbol
   */
  public Optional<GeneReport> findGeneReport(DataSource source, String geneSymbol) {
    if (m_geneReports.containsKey(source)) {
      return Optional.ofNullable(m_geneReports.get(source).get(geneSymbol));
    }
    return Optional.empty();
  }


  /**
   * Writes out {@link Phenotyper} data.
   *
   * @param outputPath the path to write a JSON file of data to
   * @throws IOException can occur from writing the file to the filesystem
   */
  public void write(Path outputPath) throws IOException {
    try (BufferedWriter writer = Files.newBufferedWriter(outputPath, StandardCharsets.UTF_8)) {
      writer.write(DataSerializer.GSON.toJson(this));
      sf_logger.info("Writing Phenotyper JSON to " + outputPath);
    }
  }

  /**
   * Read in {@link Phenotyper} data.
   */
  public static Phenotyper read(Path filePath) throws IOException {
    Preconditions.checkNotNull(filePath);
    Preconditions.checkArgument(filePath.toFile().exists());
    Preconditions.checkArgument(filePath.toFile().isFile());
    try (BufferedReader reader = Files.newBufferedReader(filePath)) {
      return DataSerializer.GSON.fromJson(reader, Phenotyper.class);
    }
  }


  private Set<String> listUnspecifiedGenes(DataSource source) {
    if (source == DataSource.UNKNOWN) {
      return Collections.emptySet();
    }
    try {
      DrugCollection cpicDrugs = new DrugCollection();
      PgkbGuidelineCollection dpwgDrugs = new PgkbGuidelineCollection();

      Set<String> unspecifiedGenes = new HashSet<>();
      if (source == DataSource.CPIC) {
        unspecifiedGenes.addAll(cpicDrugs.getAllReportableGenes());
      } else {
        unspecifiedGenes.addAll(dpwgDrugs.getGenes());
      }
      m_geneReports.get(source).values().stream()
          .map(GeneReport::getGene)
          .forEach(unspecifiedGenes::remove);
      return unspecifiedGenes;
    } catch (IOException ex) {
      throw new RuntimeException("Error reading drug data", ex);
    }
  }
}
