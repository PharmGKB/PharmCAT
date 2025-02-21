package org.pharmgkb.pharmcat.reporter.format;

import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import org.apache.commons.lang3.StringUtils;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.BaseConfig;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.PipelineWrapper;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;

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
    assertTrue(StringUtils.isNotBlank(dpydRow[13]));
    assertTrue(StringUtils.isNotBlank(dpydRow[14]));
    assertTrue(StringUtils.isNotBlank(dpydRow[15]));

    String[] cyp2c19Row = geneMap.get("CYP2C19").get(0).split("\t");
    // will have phenotype and activity score
    assertTrue(StringUtils.isNotBlank(cyp2c19Row[2]));
    assertTrue(StringUtils.isBlank(cyp2c19Row[3]));
    // but no recommendation
    assertTrue(StringUtils.isBlank(cyp2c19Row[13]));
    assertTrue(StringUtils.isBlank(cyp2c19Row[14]));
    assertTrue(StringUtils.isBlank(cyp2c19Row[15]));
  }

  @Test
  void withDebugInfo(TestInfo testInfo) throws Exception {

    try {
      System.setProperty("PHARMCAT_REPORTER_DEBUG", "true");
      PipelineWrapper testWrapper = new PipelineWrapper(testInfo, true)
          .saveIntermediateFiles();
      testWrapper.getVcfBuilder()
          .allowUnknownAllele()
          // uncallable
          .variation("TPMT", "rs1256618794", "A", "A")
          // unknown allele
          .variation("VKORC1", "rs9923231", "C", "G")
          // unknown allele, but treated as reference
          .variation("NUDT15", "rs186364861", "A", "T")
          .reference("CYP4F2")
          .reference("CYP2B6")
          .reference("CACNA1S")
          .variation("RYR1", "rs193922749", "A", "C") // c.152C>A
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

      String[] genes = {
          "CACNA1S", "CYP2B6", "CYP2C19", "CYP4F2", "DPYD", "NUDT15", "RYR1", "TPMT", "VKORC1"
      };
      List<String> referenceCalls = List.of("CYP4F2", "CYP2B6", "CACNA1S");
      List<String> noCalls = List.of("TPMT", "VKORC1");
      List<String> undocumented = List.of("NUDT15", "VKORC1");
      List<String> missingPositions = List.of("CYP2C19");

      String basename = BaseConfig.getBaseFilename(Objects.requireNonNull(vcfFile).getFileName());
      Path normalFile = vcfFile.getParent().resolve(basename + BaseConfig.REPORTER_SUFFIX + ".tsv");
      String normalTsv = Files.readString(normalFile);

      System.out.println("normal");
      System.out.println(normalTsv);
      String[] lines = normalTsv.split("\n");
      assertEquals(10, lines.length);
      Map<String, List<String>> geneMap = parseTsv(lines, 18, genes);
      assertEquals(1, geneMap.get("CYP2C19").size());
      assertEquals(1, geneMap.get("DPYD").size());

      for (String gene : genes) {
        String[] data = geneMap.get(gene).get(0).split("\t");
        if (referenceCalls.contains(gene)) {
          // check reference
          assertTrue(StringUtils.isBlank(data[12]));
          continue;
        }
        assertTrue(StringUtils.isNotBlank(data[12]));

        if (noCalls.contains(gene)) {
          // check no calls
          assertEquals("no call", data[1], "Expecting no call for " + gene);
        }

        if (missingPositions.contains(gene)) {
          assertTrue(StringUtils.isNotBlank(data[13]));
        } else {
          assertTrue(StringUtils.isBlank(data[13]));
        }

        if (undocumented.contains(gene)) {
          assertTrue(StringUtils.isNotBlank(data[14]));
          if (NamedAlleleMatcher.TREAT_UNDOCUMENTED_VARIATIONS_AS_REFERENCE.contains(gene)) {
            assertTrue(data[14].contains("treat as reference"));
          } else {
            assertFalse(data[14].contains("treat as reference"));
          }
        }

      }


    } finally {
      System.setProperty("PHARMCAT_REPORTER_DEBUG", "");
    }
  }


  @Test
  void withSampleData(TestInfo testInfo) throws Exception {
    Path sampleDataFile = TestUtils.createTempFile(testInfo, "tsv");
    try (PrintWriter writer = new PrintWriter(Files.newBufferedWriter(sampleDataFile))) {
      writer.println("PharmCAT\tTown\tStanford");
      writer.println("PharmCAT\tState\tCA");
    }

    PipelineWrapper testWrapper = new PipelineWrapper(testInfo, false)
        .saveIntermediateFiles();
    Path vcfFile = testWrapper.execute(PathUtils.getPathToResource("org/pharmgkb/pharmcat/reference.vcf"),
        null, sampleDataFile, false);

    String basename = BaseConfig.getBaseFilename(Objects.requireNonNull(vcfFile).getFileName());
    Path normalFile = testWrapper.getOutputDir().resolve(basename + BaseConfig.REPORTER_SUFFIX + ".tsv");
    String normalTsv = Files.readString(normalFile);
    String[] lines = normalTsv.split("\n");

    Map<String, List<String>> geneMap = parseTsv(lines, 18);
    Env env = new Env();
    SortedSet<String> allGenes = new TreeSet<>(env.getDefinitionReader().getGenes());
    allGenes.remove("CYP2D6");
    // TODO(markwoon): ignore NAT2 until it's fully integrated
    allGenes.remove("NAT2");
    assertEquals(allGenes, geneMap.keySet());

    System.out.println(geneMap.get("TPMT"));
    assertTrue(geneMap.get("TPMT").get(0).contains("Stanford"));
  }


  private SortedMap<String, List<String>> parseTsv(String[] lines, int maxColumns, String... genes) {
    // test normal
    SortedMap<String, List<String>> geneMap = new TreeMap<>();
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
