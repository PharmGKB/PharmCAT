package org.pharmgkb.pharmcat.phenotype;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import com.google.common.base.Preconditions;
import com.google.gson.annotations.Expose;
import com.google.gson.annotations.SerializedName;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Metadata;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.MessageAnnotation;
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

  @SerializedName("matcherMetadata")
  @Expose
  private Metadata m_matcherMetadata;
  @Expose
  @SerializedName("geneReports")
  private final SortedMap<String, GeneReport> m_geneReports = new TreeMap<>();
  @Expose
  @SerializedName("unannotatedGeneCalls")
  private SortedSet<GeneReport> m_unannotatedGeneCalls = new TreeSet<>();


  /**
   * Public constructor. This needs {@link GeneCall} objects from the {@link NamedAlleleMatcher} and {@link OutsideCall}
   * objects coming from other allele calling sources. This relies on reading definition files as well.
   *
   * @param matcherMetadata metadata for the named allele matcher used for {@code geneCalls};
   * can be null if all outside calls
   * @param geneCalls a List of {@link GeneCall} objects
   * @param outsideCalls a List of {@link OutsideCall} objects
   * @param variantWarnings map of VCF warnings, keyed to chromosomal position
   */
  public Phenotyper(Env env, @Nullable Metadata matcherMetadata, List<GeneCall> geneCalls, Set<OutsideCall> outsideCalls,
      @Nullable Map<String, Collection<String>> variantWarnings) {
    List<String> unusedGenes = initialize(geneCalls, outsideCalls, env, DataSource.CPIC, variantWarnings);
    unusedGenes.retainAll(initialize(geneCalls, outsideCalls, env, DataSource.DPWG, variantWarnings));

    if (!unusedGenes.isEmpty()) {
      for (String gene : unusedGenes) {
        GeneCall geneCall = geneCalls.stream()
            .filter(gc -> gc.getGene().equals(gene))
            .findFirst()
            .orElseThrow(() -> new IllegalStateException("Cannot find gene call for " + gene));
        GeneReport geneReport = new GeneReport(geneCall, env);
        if (!geneReport.isNoData()) {
          m_unannotatedGeneCalls.add(geneReport);
        }
      }
    }

    m_matcherMetadata = matcherMetadata;
  }


  private List<String> initialize(List<GeneCall> geneCalls, Set<OutsideCall> outsideCalls, Env env, DataSource source,
      @Nullable Map<String, Collection<String>> variantWarnings) {

    List<String> unusedGeneCalls = new ArrayList<>();
    // matcher calls
    for (GeneCall geneCall : geneCalls) {
      if (!env.hasGene(source, geneCall.getGene())) {
        unusedGeneCalls.add(geneCall.getGene());
        continue;
      }
      GeneReport geneReport = new GeneReport(geneCall, env);
      m_geneReports.put(geneReport.getGene(), geneReport);
    }

    //  outside calls
    for (OutsideCall outsideCall : outsideCalls) {
      GeneReport geneReport = m_geneReports.get(outsideCall.getGene());
      MessageAnnotation msgAnnotation = null;
      if (geneReport != null) {
        if (geneReport.getCallSource() != CallSource.OUTSIDE) {
          // outside call trumps the matcher's result
          // warn the user of the conflict
          String matcherCall = geneReport.getSourceDiplotypes().stream()
              .sorted()
              .map(Diplotype::getLabel)
              .collect(Collectors.joining("; "));
          msgAnnotation = new MessageAnnotation(MessageAnnotation.TYPE_NOTE, "prefer-sample-data",
              "PharmCAT would have called " + matcherCall + " based on the VCF but it has been ignored in " +
                  "favor of the outside call. If you want to use the PharmCAT call for this gene then remove the " +
                  "gene from the outside call data."
          );

        } else {
          // add alternate outside call
          geneReport.addOutsideCall(outsideCall, env);
          continue;
        }
      }

      geneReport = new GeneReport(outsideCall, env);
      if (msgAnnotation != null) {
        geneReport.addMessage(msgAnnotation);
      }
      m_geneReports.put(geneReport.getGene(), geneReport);
    }

    Set<String> unspecifiedGenes = listUnspecifiedGenes(env, source);
    // all other genes
    for (String geneSymbol : unspecifiedGenes) {
      m_geneReports.put(geneSymbol, GeneReport.unspecifiedGeneReport(geneSymbol, env, source));
    }

    // add VCF warnings
    m_geneReports.values().forEach(geneReport -> geneReport.addVariantWarningMessages(variantWarnings));

    return unusedGeneCalls;
  }


  public SortedMap<String, GeneReport> getGeneReports() {
    return m_geneReports;
  }

  /**
   * Find a {@link GeneReport} based on the gene symbol
   */
  public Optional<GeneReport> findGeneReport(String geneSymbol) {
    return Optional.ofNullable(m_geneReports.get(geneSymbol));
  }


  public SortedSet<GeneReport> getUnannotatedGeneCalls() {
    return m_unannotatedGeneCalls;
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
    Preconditions.checkArgument(Files.isRegularFile(filePath), "Not a file: " + filePath);
    try (BufferedReader reader = Files.newBufferedReader(filePath)) {
      return DataSerializer.GSON.fromJson(reader, Phenotyper.class);
    }
  }


  private Set<String> listUnspecifiedGenes(Env env, DataSource source) {
    if (source == DataSource.UNKNOWN) {
      return Collections.emptySet();
    }
    Set<String> unspecifiedGenes = new HashSet<>(env.getDrugs().getGenesUsedInSource(source));
    m_geneReports.values().stream()
        .map(GeneReport::getGene)
        .forEach(unspecifiedGenes::remove);
    return unspecifiedGenes;
  }


  public Metadata getMatcherMetadata() {
    return m_matcherMetadata;
  }
}
