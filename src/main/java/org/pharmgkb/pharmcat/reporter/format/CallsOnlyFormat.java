package org.pharmgkb.pharmcat.reporter.format;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.OpenOption;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.BaseConfig;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.haplotype.ResultSerializer;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.pharmgkb.pharmcat.Constants.isLowestFunctionGene;


/**
 * Generates a .tsv file that only contains calls from the Named Allele Matcher.
 *
 * @author Mark Woon
 */
public class CallsOnlyFormat extends AbstractFormat {
  private boolean m_singleFileMode;


  public CallsOnlyFormat(Path outputPath, Env env) {
    super(outputPath, env);
  }

  /**
   * Sets whether all results should be appended to a single file.
   */
  public CallsOnlyFormat singleFileMode() {
    m_singleFileMode = true;
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
      GeneReport cpicReport = reportContext.getGeneReport(DataSource.CPIC, gene);
      GeneReport dpwgReport = reportContext.getGeneReport(DataSource.DPWG, gene);
      if ((cpicReport == null || cpicReport.isNoData()) && (dpwgReport == null || dpwgReport.isNoData())) {
        continue;
      }

      GeneReport primary = cpicReport == null ? dpwgReport : cpicReport;
      calledGenes.put(gene, primary);
    }


    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(getOutputPath(), StandardCharsets.UTF_8, options))) {
      if (printHeaders) {
        writer.println("Gene\tSource Diplotype\tPhenotype\tActivity Score" +
            "\tHaplotype 1\tHaplotype 1 Function\tHaplotype 1 Activity Value" +
            "\tHaplotype 2\tHaplotype 2 Function\tHaplotype 2 Activity Value" +
            "\tOutside Call\tMatch Score\tMissing positions?\t" +
            "Recommendation Lookup Diplotype\tRecommendation Lookup Phenotype\tRecommendation Lookup Activity Score");
      }

      for (String gene : calledGenes.keySet()) {
        GeneReport report = calledGenes.get(gene);
        if (!report.isCalled()) {
          writeNoCall(writer, report);
          continue;
        }

        // only have component haplotypes for lowest function genes when diplotypes are true diplotypes
        // (vs. individual haplotypes)
        boolean lowestFunctionSingles = isLowestFunctionGene(gene) && report.getMatcherComponentHaplotypes().isEmpty();

        if (report.getSourceDiplotypes().size() > 1 || lowestFunctionSingles) {
          writeCollapsedDiplotypes(writer, report, lowestFunctionSingles, true);
        } else {
          for (Diplotype dip : report.getSourceDiplotypes()) {
            writeDiplotype(writer, report, dip, true);
          }
        }
      }

      for (GeneReport report : reportContext.getUnannotatedGeneCalls()) {
        if (report.getSourceDiplotypes().size() > 1) {
          writeCollapsedDiplotypes(writer, report, false, false);
        } else {
          for (Diplotype dip : report.getSourceDiplotypes()) {
            writeDiplotype(writer, report, dip, false);
          }
        }
      }
    }
  }

  private void writeNoCall(PrintWriter writer, GeneReport report) {
    writer.print(report.getGene());
    writer.print("\tno call\t\t" +
        "\t\t\t" +
        "\t\t\t");
    writeCommon(writer, report, null, false);
    writer.println();
  }


  private boolean isIgnorableValue(String text) {
    return StringUtils.isBlank(text) || text.equals(TextConstants.NA) || text.equals(TextConstants.NO_RESULT);
  }

  private String generateStandardizedValue(String text) {
    return isIgnorableValue(text) ? " " : text;
  }

  private String generatePhenotypeValue(List<String> phenotypes) {
    return phenotypes.stream()
        .map(this::generateStandardizedValue)
        .collect(Collectors.joining(", "));
  }


  private void writeCollapsedDiplotypes(PrintWriter writer, GeneReport report, boolean lowestFunctionSingles,
      boolean showRecommendationDiplotype) {

    boolean hasPhenotypes = report.getSourceDiplotypes().stream()
        .anyMatch(d -> !d.getPhenotypes().isEmpty() && !isIgnorableValue(d.getPhenotypes().get(0)));
    boolean hasActivityScores = report.getSourceDiplotypes().stream()
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
    };

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

    writeCommon(writer, report, matchScores.toString(), showRecommendationDiplotype);
    writer.println();
  }

  private void writeDiplotype(PrintWriter writer, GeneReport report, Diplotype dip,
      boolean showRecommendationDiplotype) {
    writer.print(report.getGene());
    writer.print("\t");
    // diplotype
    writer.print(buildDiplotypeName(dip, report));
    writer.print("\t");
    // phenotype
    writer.print(generatePhenotypeValue(dip.getPhenotypes()));
    writer.print("\t");
    // activity score
    if (dip.getActivityScore() != null) {
      writer.print(generateStandardizedValue(dip.getActivityScore()));
    }
    writer.print("\t");
    // haplotype 1
    if (dip.getAllele1() != null) {
      writer.print(dip.getAllele1().getName());
      writer.print("\t");
      if (dip.getAllele1().getFunction() != null) {
        writer.print(dip.getAllele1().getFunction());
      }
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
      writer.print(dip.getAllele2());
      writer.print("\t");
      if (dip.getAllele2().getFunction() != null) {
        writer.print(dip.getAllele2().getFunction());
      }
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

    writeCommon(writer, report, Integer.toString(dip.getMatchScore()), showRecommendationDiplotype);
    writer.println();
  }


  private void writeCommon(PrintWriter writer, GeneReport report, String matchScore,
      boolean showRecommendationDiplotype) {
    // outside call
    writer.print(report.isOutsideCall());
    writer.print("\t");
    writer.print(matchScore);
    writer.print("\t");
    // missing positions
    writer.print(report.isMissingVariants());
    writer.print("\t");
    // recommendation lookup fields
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
      writer.print(generateStandardizedValue(recDip.getActivityScore()));
    } else {
      writer.print("\t\t");
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


  public static void main(String[] args) {

    if (args == null || args.length != 2) {
      throw new IllegalArgumentException("Please specify an input and output directory");
    }
    Path inDir = Paths.get(args[0]);
    if (!Files.isDirectory(inDir)) {
      throw new IllegalArgumentException("Not a valid directory: " + inDir);
    }
    Path outDir = Paths.get(args[1]);
    if (!Files.isDirectory(outDir) && Files.exists(outDir)) {
      throw new IllegalArgumentException("Not a valid directory: " + outDir);
    }

    System.out.println("Reading from " + inDir);
    System.out.println("Writing to " + outDir);

    try {
      if (!Files.exists(outDir)) {
        Files.createDirectories(outDir);
      }

      Env env = new Env();
      readDir(env, inDir, outDir, 0);

    } catch (Exception ex) {
      //noinspection CallToPrintStackTrace
      ex.printStackTrace();
    }
  }

  private static int readDir(Env env, Path inDir, Path outDir, int index) throws IOException {

    List<Path> matchFiles = new ArrayList<>();
    List<Path> dirs = new ArrayList<>();

    try (Stream<Path> files = Files.list(inDir)) {
      files.forEach(f -> {
        if (Files.isDirectory(f)) {
          dirs.add(f);
        } else if (Files.isRegularFile(f)) {
          if (f.toString().endsWith(".match.json")) {
            matchFiles.add(f);
          }
        }
      });
    }

    System.out.println("Found " + matchFiles.size() + " in " + inDir);
    for (Path d : dirs) {
      index = readDir(env, d, outDir, index);
    }
    for (Path mFile : matchFiles) {
      index += 1;
      //String basename = BaseConfig.getBaseFilename(mFile);
      String basename = "complete";
      if (index % 1000 == 0) {
        System.out.println(index);
      }
      org.pharmgkb.pharmcat.haplotype.model.Result matcherResult = new ResultSerializer()
          .fromJson(mFile);
      Phenotyper phenotyper = new Phenotyper(env, matcherResult.getMetadata(), matcherResult.getGeneCalls(),
          new TreeSet<>(), new HashMap<>());

      ReportContext reportContext = new ReportContext(env, phenotyper, basename);

      Path outFile = outDir.resolve(basename + BaseConfig.REPORTER_SUFFIX + ".tsv");
      new CallsOnlyFormat(outFile, env)
          .singleFileMode()
          .write(reportContext);
    }
    return index;
  }
}
