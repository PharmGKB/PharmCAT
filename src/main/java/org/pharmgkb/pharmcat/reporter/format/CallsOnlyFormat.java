package org.pharmgkb.pharmcat.reporter.format;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.stream.Collectors;
import org.apache.commons.lang3.StringUtils;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.VariantReport;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;
import org.pharmgkb.pharmcat.util.CliUtils;

import static org.pharmgkb.pharmcat.Constants.isActivityScoreGene;
import static org.pharmgkb.pharmcat.Constants.isLowestFunctionGene;


/**
 * Generates a .tsv file that only contains calls from the Named Allele Matcher.
 *
 * @author Mark Woon
 */
public class CallsOnlyFormat extends AbstractFormat {
  public static final String ENV_DEBUG_KEY = "PHARMCAT_REPORTER_DEBUG";
  public static final String NO_CALL_TAG = "no call";
  public static final String TREAT_AS_REFERENCE_TAG = "(treat as reference)";
  public static final String HEADER_SAMPLE_ID = "Sample ID";
  public static final String HEADER_VARIANTS = "Variants";
  public static final String HEADER_UNDOCUMENTED_VARIANTS = "Undocumented variants";
  private boolean m_singleFileMode;
  // this is only used in single file mode
  private boolean m_showSampleId = true;
  private boolean m_showVariants;
  private boolean m_showMissingVariants;
  private boolean m_showUndocumentedVariants;
  private final boolean m_debug;


  public CallsOnlyFormat(Path outputPath, Env env) {
    super(outputPath, env);
    m_debug = Boolean.parseBoolean(System.getenv(ENV_DEBUG_KEY)) ||
        Boolean.parseBoolean(System.getProperty(ENV_DEBUG_KEY));
    if (m_debug) {
      m_showVariants = true;
      m_showMissingVariants = true;
      m_showUndocumentedVariants = true;
    }
  }


  /**
   * Sets whether results should include variants.
   */
  public CallsOnlyFormat showVariants() {
    m_showVariants = true;
    return this;
  }


  /**
   * Sets whether results should specify missing variants (defaults to yes/no).
   */
  public CallsOnlyFormat showMissingVariants() {
    m_showMissingVariants = true;
    return this;
  }

  /**
   * Sets whether results should specify undocumented variants.
   */
  public CallsOnlyFormat showUndocumentedVariants() {
    m_showUndocumentedVariants = true;
    return this;
  }

  /**
   * Sets whether all results should be appended to a single file.
   */
  public CallsOnlyFormat singleFileMode() {
    m_singleFileMode = true;
    return this;
  }

  /**
   * Sets whether the Sample ID column should be written.
   * This only matters in single-file mode.
   */
  public CallsOnlyFormat hideSampleId() {
    m_showSampleId = false;
    return this;
  }


  @Override
  public void write(ReportContext reportContext) throws IOException {

    boolean printHeaders = true;
    OpenOption[] options = new OpenOption[] {
        StandardOpenOption.CREATE,
        StandardOpenOption.WRITE,
        StandardOpenOption.TRUNCATE_EXISTING,
    };
    if (m_singleFileMode && Files.exists(getOutputPath())) {
      printHeaders = false;
      options = new OpenOption[] {
          StandardOpenOption.CREATE,
          StandardOpenOption.WRITE,
          StandardOpenOption.APPEND,
      };
    }

    Map<String, GeneReport> calledGenes = new HashMap<>();
    for (String gene : getEnv().getDefinitionReader().getGenes()) {
      GeneReport geneReport = reportContext.getGeneReport(gene);
      if (geneReport == null || geneReport.isNoData()) {
        continue;
      }

      calledGenes.put(gene, geneReport);
    }

    String sampleId = null;
    SortedMap<String, String> sampleProps = null;
    // older versions of reports do not have metadata
    //noinspection ConstantValue
    if (reportContext.getMatcherMetadata() != null) {
      sampleId = reportContext.getMatcherMetadata().getSampleId();
      if (reportContext.getMatcherMetadata().getSampleProps() != null &&
          !reportContext.getMatcherMetadata().getSampleProps().isEmpty()) {
        sampleProps = new TreeMap<>(reportContext.getMatcherMetadata().getSampleProps());
      }
    }
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(getOutputPath(), StandardCharsets.UTF_8, options))) {
      if (printHeaders) {
        writer.print("PharmCAT " + CliUtils.getVersion());
        writer.println();
        if (m_singleFileMode && m_showSampleId) {
          writer.print(HEADER_SAMPLE_ID + "\t");
        }
        writer.print("Gene\tSource Diplotype\tPhenotype\tActivity Score" +
            "\tHaplotype 1\tHaplotype 1 Function\tHaplotype 1 Activity Value" +
            "\tHaplotype 2\tHaplotype 2 Function\tHaplotype 2 Activity Value" +
            "\tOutside Call\tMatch Score\t");
        if (m_showVariants) {
          writer.print(HEADER_VARIANTS + "\t");
        }
        writer.print("Missing positions\t");
        if (m_showUndocumentedVariants) {
          writer.print(HEADER_UNDOCUMENTED_VARIANTS + "\t");
        }
        writer.print("Recommendation Lookup Diplotype\tRecommendation Lookup Phenotype\t" +
            "Recommendation Lookup Activity Score");
        if (sampleProps != null) {
          for (String key : sampleProps.keySet()) {
            writer.print("\t");
            writer.print(key);
          }
        }
        writer.println();
      }

      for (String gene : calledGenes.keySet()) {
        GeneReport report = calledGenes.get(gene);
        if (!report.isCalled()) {
          writeNoCall(writer, sampleId, sampleProps, report);
          continue;
        }


        // only have component haplotypes for lowest function genes when diplotypes are true diplotypes
        // (vs. individual haplotypes)
        boolean lowestFunctionSingles = isLowestFunctionGene(gene) && report.getMatcherComponentHaplotypes().isEmpty();

        if (report.getSourceDiplotypes().size() > 1 || lowestFunctionSingles) {
          writeCollapsedDiplotypes(writer, sampleId, sampleProps, report, lowestFunctionSingles, true);
        } else {
          for (Diplotype dip : report.getSourceDiplotypes()) {
            writeDiplotype(writer, sampleId, sampleProps, report, dip, true);
          }
        }
      }

      for (GeneReport report : reportContext.getUnannotatedGeneCalls()) {
        if (report.getSourceDiplotypes().size() > 1) {
          writeCollapsedDiplotypes(writer, sampleId, sampleProps, report, false, false);
        } else {
          for (Diplotype dip : report.getSourceDiplotypes()) {
            writeDiplotype(writer, sampleId, sampleProps, report, dip, false);
          }
        }
      }
    }
  }

  private void writeNoCall(PrintWriter writer, @Nullable String sampleId, @Nullable Map<String, String> sampleProps,
      GeneReport report) {
    if (m_singleFileMode && m_showSampleId) {
      if (sampleId != null) {
        writer.print(sampleId);
      }
      writer.print("\t");
    }
    writer.print(report.getGene());
    writer.print("\t" + NO_CALL_TAG + "\t\t\t" +
        "\t\t\t" +
        "\t\t\t");
    writeCommon(writer, sampleProps, report, null, false);
    writer.println();
  }


  private boolean isIgnorableValue(@Nullable String text) {
    return StringUtils.isBlank(text) || text.equals(TextConstants.NA) || text.equals(TextConstants.NO_RESULT);
  }

  private String generateStandardizedValue(@Nullable String text) {
    return isIgnorableValue(text) ? " " : text;
  }

  private String generatePhenotypeValue(List<String> phenotypes) {
    return phenotypes.stream()
        .map(this::generateStandardizedValue)
        .collect(Collectors.joining(", "));
  }


  private void writeCollapsedDiplotypes(PrintWriter writer,
      @Nullable String sampleId, @Nullable Map<String, String> sampleProps,
      GeneReport report, boolean lowestFunctionSingles, boolean showRecommendationDiplotype) {

    boolean hasPhenotypes = !lowestFunctionSingles && report.getSourceDiplotypes().stream()
        .anyMatch(d -> !d.getPhenotypes().isEmpty() && !isIgnorableValue(d.getPhenotypes().get(0)));
    boolean hasActivityScores = isActivityScoreGene(report.getGene()) &&
        report.getSourceDiplotypes().stream()
            .anyMatch(d -> !isIgnorableValue(d.getActivityScore()));

    StringBuilder diplotypes = new StringBuilder();
    StringBuilder matchScores = new StringBuilder();
    StringBuilder phenotypes = new StringBuilder();
    StringBuilder activityScores = new StringBuilder();
    for (Diplotype dip : report.getSourceDiplotypes()) {
      if (!diplotypes.isEmpty()) {
        diplotypes.append(lowestFunctionSingles ? " AND " : " OR ");
      }
      diplotypes.append(buildDiplotypeName(dip, report));

      if (hasPhenotypes) {
        if (!phenotypes.isEmpty()) {
          phenotypes.append(" / ");
        }
        phenotypes.append(generatePhenotypeValue(dip.getPhenotypes()));
      }

      if (hasActivityScores) {
        if (!activityScores.isEmpty()) {
          activityScores.append(" / ");
        }
        activityScores.append(generateStandardizedValue(dip.getActivityScore()));
      }

      if (!lowestFunctionSingles) {
        if (!matchScores.isEmpty()) {
          matchScores.append(" / ");
        }
        matchScores.append(dip.getMatchScore());
      }
    }

    if (m_singleFileMode && m_showSampleId) {
      if (sampleId != null) {
        writer.print(sampleId);
      }
      writer.print("\t");
    }
    writer.print(report.getGene());
    writer.print("\t");
    writer.print(diplotypes);
    writer.print("\t");
    writer.print(phenotypes);
    writer.print("\t");
    writer.print(activityScores);
    writer.print("\t" +
        "\t\t\t" +
        "\t\t\t");

    writeCommon(writer, sampleProps, report, matchScores.toString(), showRecommendationDiplotype);
    writer.println();
  }

  private void writeDiplotype(PrintWriter writer, @Nullable String sampleId, @Nullable Map<String, String> sampleProps,
      GeneReport report, Diplotype dip, boolean showRecommendationDiplotype) {

    if (m_singleFileMode && m_showSampleId) {
      if (sampleId != null) {
        writer.print(sampleId);
      }
      writer.print("\t");
    }
    writer.print(report.getGene());
    writer.print("\t");
    // diplotype
    writer.print(buildDiplotypeName(dip, report));
    writer.print("\t");
    // phenotype
    writer.print(generatePhenotypeValue(dip.getPhenotypes()));
    writer.print("\t");
    // activity score
    if (isActivityScoreGene(report.getGene()) && dip.getActivityScore() != null) {
      writer.print(generateStandardizedValue(dip.getActivityScore()));
    }
    writer.print("\t");
    // haplotype 1
    if (dip.getAllele1() != null) {
      writer.print(dip.getAllele1().getName());
      writer.print("\t");
      writer.print(dip.getAllele1().getFunction());
      writer.print("\t");
      if (dip.getAllele1().getActivityValue() != null &&
          !dip.getAllele1().getActivityValue().equals(TextConstants.NA)) {
        writer.print(dip.getAllele1().getActivityValue());
      }
    } else {
      writer.print("\t");
      writer.print("\t");
    }
    writer.print("\t");
    // haplotype 2
    if (dip.getAllele2() != null) {
      writer.print(dip.getAllele2().getName());
      writer.print("\t");
      writer.print(dip.getAllele2().getFunction());
      writer.print("\t");
      if (dip.getAllele2().getActivityValue() != null &&
          !dip.getAllele2().getActivityValue().equals(TextConstants.NA)) {
        writer.print(dip.getAllele2().getActivityValue());
      }
    } else {
      writer.print("\t");
      writer.print("\t");
    }
    writer.print("\t");

    writeCommon(writer, sampleProps, report, Integer.toString(dip.getMatchScore()), showRecommendationDiplotype);
    writer.println();
  }


  private void writeCommon(PrintWriter writer,  @Nullable Map<String, String> sampleProps, GeneReport report,
      @Nullable String matchScore, boolean showRecommendationDiplotype) {
    // outside call
    writer.print(report.isOutsideCall() ? "yes" : "no");
    writer.print("\t");
    writer.print(matchScore == null ? "" : matchScore);
    writer.print("\t");
    // variants
    if (m_showVariants) {
      StringBuilder varBuilder = new StringBuilder();
      for (VariantReport vr : report.getVariantReports()) {
        if (!vr.isMissing() && vr.isNonReference()) {
            if (!varBuilder.isEmpty()) {
              varBuilder.append(", ");
            }
            varBuilder.append(vr.getPosition())
                .append(":")
                .append(vr.getCall());
          }
        }
      writer.print(varBuilder);
      writer.print("\t");
    }
    // missing positions
    if (m_showMissingVariants) {
      StringBuilder missingBuilder = new StringBuilder();
      for (VariantReport vr : report.getVariantReports()) {
        if (vr.isMissing()) {
          if (m_showMissingVariants) {
            if (!missingBuilder.isEmpty()) {
              missingBuilder.append(", ");
            }
            missingBuilder.append(vr.getPosition());
          }
        }
      }
      writer.print(missingBuilder);
    } else {
      writer.print(report.isMissingVariants() ? "yes" : "no");
    }
    writer.print("\t");
    if (m_showUndocumentedVariants) {
      if (report.isHasUndocumentedVariations()) {
        writer.print(report.getVariantReports().stream()
            .filter(VariantReport::isHasUndocumentedVariations)
            .map(vr -> Long.toString(vr.getPosition()))
            .collect(Collectors.joining(", ")));
        if (report.isTreatUndocumentedVariationsAsReference()) {
          writer.print(" " + TREAT_AS_REFERENCE_TAG);
        }
      }
      writer.print("\t");
    }
    // recommendation lookup fields
    //noinspection ConstantValue
    if (showRecommendationDiplotype && report.getRecommendationDiplotypes() != null &&
        report.getRecommendationDiplotypes().size() == 1) {
      Diplotype recDip = report.getRecommendationDiplotypes().first();
      // recommendation lookup diplotype
      writer.print(buildDiplotypeName(recDip, report));
      writer.print("\t");
      // recommendation lookup phenotype
      writer.print(generatePhenotypeValue(recDip.getPhenotypes()));
      writer.print("\t");
      // recommendation lookup activity score
      if (isActivityScoreGene(report.getGene())) {
        writer.print(generateStandardizedValue(Objects.requireNonNull(recDip.getActivityScore())));
      }
    } else {
      writer.print("\t\t");
    }
    if (sampleProps != null) {
      for (String key : sampleProps.keySet()) {
        writer.print("\t");
        writer.print(sampleProps.get(key));
      }
    }
  }

  private String buildDiplotypeName(Diplotype dip, GeneReport geneReport) {
    StringBuilder builder = new StringBuilder();
    if (dip.getAllele1() != null) {
      builder.append(dip.getAllele1().getName());
      if (dip.getAllele2() != null) {
        builder.append("/")
            .append(dip.getAllele2().getName());
      } else if (geneReport.getMatcherHomozygousComponentHaplotypes().contains(dip.getAllele1().getName())) {
        builder.append("/")
            .append(dip.getAllele1().getName());
      }
    }
    return builder.toString();
  }
}
