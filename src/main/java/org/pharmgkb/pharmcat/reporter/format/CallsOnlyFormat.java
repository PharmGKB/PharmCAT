package org.pharmgkb.pharmcat.reporter.format;

import java.io.IOException;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import org.pharmgkb.pharmcat.BaseConfig;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.Phenotyper;
import org.pharmgkb.pharmcat.reporter.ReportContext;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;


/**
 * Generates a .tsv file that only contains calls from the Named Allele Matcher.
 *
 * @author Mark Woon
 */
public class CallsOnlyFormat extends AbstractFormat {
  private boolean m_showMatchScores;


  public CallsOnlyFormat(Path outputPath, Env env) {
    super(outputPath, env);
  }

  /**
   * Sets whether match scores should be exported.
   * This is only necessary if not calling with top-candidates-only.
   */
  public CallsOnlyFormat showMatchScores() {
    m_showMatchScores = true;
    return this;
  }


  @Override
  public void write(ReportContext reportContext) throws IOException {
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(getOutputPath(), StandardCharsets.UTF_8))) {
      writer.print("Gene\tDiplotype\tPhenotype\tActivity Score" +
          "\tHaplotype 1\tHaplotype 1 Function\tHaplotype 1 Activity Value" +
          "\tHaplotype 2\tHaplotype 2 Function\tHaplotype 2 Activity Value" +
          "\tOutside Call\t");
      if (m_showMatchScores) {
        writer.print("Match Score\t");
      }
      writer.println("Missing positions?");

      for (String gene : getEnv().getDefinitionReader().getGenes()) {
        GeneReport cpicReport = reportContext.getGeneReport(DataSource.CPIC, gene);
        GeneReport dpwgReport = reportContext.getGeneReport(DataSource.DPWG, gene);
        if ((cpicReport == null || !cpicReport.isCalled()) && (dpwgReport == null || !dpwgReport.isCalled())) {
          continue;
        }

        GeneReport primary = (cpicReport == null || !cpicReport.isCalled()) ? dpwgReport : cpicReport;
        for (Diplotype dip : primary.getSourceDiplotypes()) {
          writer.print(gene);
          writer.print("\t");
          // diplotype
          if (dip.getAllele1() != null) {
            writer.print(dip.getAllele1().getName());
            if (dip.getAllele2() != null) {
              writer.print("/");
              writer.print(dip.getAllele2().getName());
            }
          }
          writer.print("\t");
          // phenotype
          if (!dip.getPhenotypes().isEmpty()) {
            writer.print(dip.getPhenotypes().stream()
                .filter(p -> !p.equals(TextConstants.NO_RESULT))
                .collect(Collectors.joining(", ")));
          }
          writer.print("\t");
          // activity score
          if (dip.getActivityScore() != null) {
            writer.print(dip.getActivityScore());
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
          // outside call
          writer.print(primary.isOutsideCall());
          writer.print("\t");
          if (m_showMatchScores) {
            // match score
            writer.print(dip.getMatchScore());
            writer.print("\t");
          }
          // missing positions
          writer.print(primary.isMissingVariants());
          writer.println();
        }
      }
    }
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

    List<Path> phenotypeFiles = new ArrayList<>();
    List<Path> dirs = new ArrayList<>();

    try (Stream<Path> files = Files.list(inDir)) {
      files.forEach(f -> {
        if (Files.isDirectory(f)) {
          dirs.add(f);
        } else if (Files.isRegularFile(f)) {
          if (f.toString().endsWith("phenotype.json")) {
            phenotypeFiles.add(f);
          }
        }
      });
    }

    System.out.println("Found " + phenotypeFiles.size() + " in " + inDir);
    for (Path d : dirs) {
      index = readDir(env, d, outDir, index);
    }
    for (Path pFile : phenotypeFiles) {
      index += 1;
      String basename = String.format("%06d", index);
      //System.out.println(pFile + " -> " + basename);
      Phenotyper phenotyper = Phenotyper.read(pFile);
      ReportContext reportContext = new ReportContext(env, phenotyper.getGeneReports(), basename);

      Path outFile = outDir.resolve(basename + BaseConfig.REPORTER_SUFFIX + ".tsv");
      new CallsOnlyFormat(outFile, env)
          .write(reportContext);
    }
    return index;
  }
}
