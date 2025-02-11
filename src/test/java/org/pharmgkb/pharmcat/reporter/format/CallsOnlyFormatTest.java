package org.pharmgkb.pharmcat.reporter.format;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import org.apache.commons.lang3.StringUtils;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.BaseConfig;
import org.pharmgkb.pharmcat.PipelineWrapper;

import static org.junit.jupiter.api.Assertions.*;
import static org.pharmgkb.pharmcat.reporter.model.result.Haplotype.UNKNOWN;


/**
 * This is a JUnit test for {@link CallsOnlyFormat}.
 *
 * @author Mark Woon
 */
class CallsOnlyFormatTest {


  @Test
  void reference(TestInfo testInfo) throws Exception {
    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false)
        .saveIntermediateFiles()
        ;
    Path vcfFile = testWrapper.executeWithVcf(PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf"));

    String basename = BaseConfig.getBaseFilename(Objects.requireNonNull(vcfFile).getFileName());
    Path normalFile = testWrapper.getOutputDir().resolve(basename + BaseConfig.REPORTER_SUFFIX + ".tsv");
    String normalTsv = Files.readString(normalFile);
    String[] lines = normalTsv.split("\n");
    Map<String, List<String>> geneMap = parseTsv(lines, 16);
    assertTrue(geneMap.size() >= 18);
  }


  @Test
  void dpydHaplotypes(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false)
        .saveIntermediateFiles();
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs72547601", "C", "C") // c.2933A>G - no function
        .variation("DPYD", "rs67376798", "A", "T") // c.2846A>T - decreased
        .variation("DPYD", "rs60139309", "T", "C") // c.2582A>G - normal
    ;
    Path vcfFile = testWrapper.execute();

    List<String> expectedCalls = List.of("c.2582A>G", "c.2846A>T", "c.2933A>G/c.2933A>G");

    String basename = BaseConfig.getBaseFilename(Objects.requireNonNull(vcfFile).getFileName());
    Path normalFile = vcfFile.getParent().resolve(basename + BaseConfig.REPORTER_SUFFIX + ".tsv");
    String normalTsv = Files.readString(normalFile);
    String[] lines = normalTsv.split("\n");
    Map<String, List<String>> geneMap = parseTsv(lines, 16, "DPYD");
    assertEquals(1, geneMap.get("DPYD").size());
    checkTextContains(geneMap.get("DPYD").get(0), expectedCalls);
  }


  @Test
  void multipleDiplotypes(TestInfo testInfo) throws Exception {

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true)
        .saveIntermediateFiles();
    testWrapper.getVcfBuilder()
        .variation("DPYD", "rs72547601", "C", "C") // c.2933A>G - no function
        .variation("DPYD", "rs67376798", "A", "T") // c.2846A>T - decreased
        .variation("DPYD", "rs60139309", "T", "C") // c.2582A>G - normal
        .missing("CYP2C19",
            "rs55752064",
            "rs1564656981",
            "rs55640102")
        .variation("CYP2C19", "rs3758581", "G", "G")
    ;
    Path vcfFile = testWrapper.execute();

    List<String> expectedDpydCalls = List.of("c.2582A>G", "c.2846A>T", "c.2933A>G/c.2933A>G");

    String basename = BaseConfig.getBaseFilename(Objects.requireNonNull(vcfFile).getFileName());
    Path normalFile = vcfFile.getParent().resolve(basename + BaseConfig.REPORTER_SUFFIX + ".tsv");
    String normalTsv = Files.readString(normalFile);

    System.out.println("normal");
    System.out.println(normalTsv);
    String[] lines = normalTsv.split("\n");
    assertEquals(3, lines.length);
    Map<String, List<String>> geneMap = parseTsv(lines, 16, "CYP2C19", "DPYD");
    assertEquals(1, geneMap.get("CYP2C19").size());
    assertEquals(1, geneMap.get("DPYD").size());

    checkTextContains(geneMap.get("DPYD").get(0), expectedDpydCalls);
    String[] dpydRow = geneMap.get("DPYD").get(0).split("\t");
    // no phenotype and activity score, but recommendation phenotype and activity score
    assertTrue(StringUtils.isBlank(dpydRow[2]));
    assertTrue(StringUtils.isBlank(dpydRow[3]));
    System.out.println(dpydRow[14]);
    System.out.println(dpydRow[15]);
    assertTrue(StringUtils.isNotBlank(dpydRow[14]));
    assertTrue(StringUtils.isNotBlank(dpydRow[15]));
  }


  private Map<String, List<String>> parseTsv(String[] lines, int maxColumns, String... genes) {
    // test normal
    Map<String, List<String>> geneMap = new HashMap<>();
    for (String line : lines) {
      String[] data = line.split("\t");
      if (data.length != maxColumns) {
        StringBuilder builder = new StringBuilder();
        for (int x = 0; x < data.length; x += 1) {
          builder.append(x)
              .append(" - ")
              .append(data[x])
              .append("\n");
        }
        assertEquals(maxColumns, data.length, "Found columns:\n" + builder);
      }
      if (line.startsWith("Gene\t")) {
        continue;
      }
      geneMap.computeIfAbsent(data[0], k -> new ArrayList<>())
          .add(line);
    }
    if (genes != null && genes.length > 0) {
      assertEquals(genes.length, geneMap.keySet().stream()
          .map(g -> geneMap.get(g).get(0))
          .filter(n -> !n.contains(UNKNOWN))
          .count());
      for (String gene : genes) {
        assertNotNull(geneMap.get(gene), "Missing data for " + gene);
      }
    }
    return geneMap;
  }

  private void checkTextContains(String text, List<String> calls) {
    for (String call : calls) {
      assertTrue(text.contains(call));
    }
  }
}
