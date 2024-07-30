package org.pharmgkb.pharmcat;

import java.io.BufferedReader;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import com.google.common.base.Splitter;
import org.apache.commons.lang3.StringUtils;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.reporter.handlebars.ReportHelpers;
import org.pharmgkb.pharmcat.reporter.model.DataSource;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.pharmgkb.pharmcat.PipelineTest.readHtmlReport;


/**
 * JUnit test to check SHC modifications.
 *
 * @author Mark Woon
 */
class ShcTest {
  private static final Splitter sf_commaSplitter = Splitter.on(",").trimResults().omitEmptyStrings();
  private static final Splitter sf_orSplitter = Splitter.on(" or ").trimResults().omitEmptyStrings();
  private static final Pattern sf_varPattern = Pattern.compile("(rs\\d+) .*\\(((\\d)([/|])(\\d))\\)");


  @BeforeAll
  static void prepare() {
    ReportHelpers.setDebugMode(true);
    //TestUtils.setSaveTestOutput(true);
  }


  @Test
  void test(TestInfo testInfo) throws Exception {
    // source data: https://docs.google.com/spreadsheets/d/1Yhl8x_xJaNaxky_YtOUU-DZJj2VV6z95trURiIcSTrU/edit
    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/ShcTest.tsv");
    List<List<String[]>> tests = new ArrayList<>();
    try (BufferedReader reader = Files.newBufferedReader(tsvFile)) {
      String line = reader.readLine();
      int x = 0;
      List<String[]> test = null;
      while (line != null) {
        if (line.startsWith("---")) {
          line = reader.readLine();
          continue;
        }

        if (test == null) {
          if (!line.contains("_Report")) {
            line = reader.readLine();
            continue;
          }
          test = new ArrayList<>();
          tests.add(test);
        }
        if (StringUtils.isBlank(line)) {
          test = null;
          line = reader.readLine();
          continue;
        }
        String[] data = line.split("\t");
        test.add(data);

        line = reader.readLine();
        x += 1;
      }
    }

    DefinitionReader definitionReader = new DefinitionReader();
    int idx = 0;
    for (List<String[]> test : tests) {
      idx += 1;
      if (idx < 10) {
        //continue;
      }
      String name = test.get(0)[0];
      System.out.println(name);

      PipelineWrapper testWrapper = new PipelineWrapper(testInfo, name, false, false, false)
          .saveIntermediateFiles();
      TestVcfBuilder vcfBuilder = testWrapper.getVcfBuilder()
          .phased();
      Path outsideCallPath = TestUtils.createTestFile(testInfo, ".tsv");

      boolean hasOutsideCalls = false;
      try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(outsideCallPath))) {
        for (int x = 1; x < test.size(); x += 1) {
          String[] data = test.get(x);
          String condition = data[0];
          if ("no input".equalsIgnoreCase(StringUtils.stripToEmpty(condition))) {
            continue;
          }

          String gene = data[1];
          if (condition.toLowerCase().startsWith("outside call")) {
            //System.out.println(gene + " " + data[2]);
            writer.println(gene + "\t" + data[2]);
            hasOutsideCalls = true;
          } else {
            //System.out.println(gene + " " + data[2] + " - " + data[3]);
            DefinitionFile definitionFile = definitionReader.getDefinitionFile(gene);
            for (String var : sf_commaSplitter.splitToList(data[3])) {
              Matcher m = sf_varPattern.matcher(var);
              if (!m.matches()) {
                System.out.println("UNSUPPORTED FORMAT: " + gene + " " + var);
                //throw new IllegalStateException("Unsupported format: " + var);
                continue;
              }
              String rsid = m.group(1);
              VariantLocus vl = Arrays.stream(definitionFile.getVariants())
                  .filter(v -> rsid.equals(v.getRsid()))
                  .findFirst()
                  .orElseThrow(() -> new IllegalStateException("Can't find '" + rsid + "'"));

              String gt = m.group(2);
              vcfBuilder.variationAsIs(gene, rsid, gt, vl.getRef(), vl.getAlts().toArray(new String[0]));
            }

            for (String var : sf_commaSplitter.splitToList(data[4])) {
              String[] undoc = var.split(" ");
              String rsid = undoc[0];
              String v1 = undoc[1];
              String v2 = undoc[2];

              VariantLocus vl = Arrays.stream(definitionFile.getVariants())
                  .filter(v -> rsid.equals(v.getRsid()))
                  .findFirst()
                  .orElseThrow(() -> new IllegalStateException("Can't find " + rsid));

              List<String> gt = new ArrayList<>();
              String[] alts = new String[1];
              if (v1.equals(vl.getRef())) {
                gt.add("0");
              } else {
                gt.add("1");
                alts[0] = v1;
              }
              if (v2.equals(vl.getRef())) {
                gt.add("0");
              } else {
                gt.add("1");
                alts[0] = v2;
              }
              vcfBuilder.variationAsIs(gene, rsid, String.join("|", gt), vl.getRef(), alts);
            }
          }
        }
      }

      Path vcfFile;
      if (hasOutsideCalls) {
        vcfFile = testWrapper.execute(outsideCallPath);
      } else {
        vcfFile = testWrapper.execute();
      }

      for (int x = 1; x < test.size(); x += 1) {
        String[] data = test.get(x);
        String condition = data[0];
        if ("no input".equalsIgnoreCase(StringUtils.stripToEmpty(condition))) {
          continue;
        }

        String gene = data[1];
        List<String> rez = sf_orSplitter.splitToList(data[5]).stream()
            //.sorted(HaplotypeNameComparator.getComparator())
            .toList();
        if (gene.equals("CYP3A4")) {
          testWrapper.testPrintCalls(DataSource.DPWG, gene, rez);
        } else {
          testWrapper.testPrintCalls(DataSource.CPIC, gene, rez);
        }
        if (!gene.equals("DPYD") && data.length > 6) {
          String pheno = StringUtils.stripToNull(data[6]);
          if (pheno != null) {
            if (pheno.contains(" AS ")) {
              pheno = data[6].split(" ")[0];
            }
            switch (pheno) {
              case "IM" -> pheno = "Intermediate Metabolizer";
              case "NM" -> pheno = "Normal Metabolizer";
              case "PM" -> pheno = "Poor Metabolizer";
              case "UM" -> pheno = "Ultrarapid Metabolizer";
            }

            Document document = readHtmlReport(vcfFile);
            Elements gsDips = document.select(".gs-" + gene + " .gs-dip td");
            if (gsDips.isEmpty()) {
              System.out.println("WTF");
            }
            assertEquals(0, gsDips.size() % 3);
            for (int y = 0; y < gsDips.size() % 3; y += 1) {
              assertEquals(pheno, gsDips.get(((x + 1) * 3) - 1).text(), "Phenotype mismatch for " + gene);
            }
          }
        }
      }

      System.out.println();
      System.out.println();
    }
  }
}
