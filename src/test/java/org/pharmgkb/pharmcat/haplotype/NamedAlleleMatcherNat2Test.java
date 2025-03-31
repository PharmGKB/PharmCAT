package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import com.google.common.base.Splitter;
import org.apache.commons.lang3.StringUtils;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.junit.jupiter.api.Assertions.*;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.printMatches;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * Test calling for NAT2.
 *
 * @author Mark Woon
 */
public class NamedAlleleMatcherNat2Test {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("NAT2_translation.json");
  private static Env s_env = null;
  private static final Splitter sf_orSplitter = Splitter.on("or").omitEmptyStrings().trimResults();

  @BeforeAll
  static void prepare() throws IOException, ReportableException {
    s_env = new Env();
    //TestUtils.setSaveTestOutput(true);
  }

  @AfterEach
  void deleteDirectory(TestInfo testInfo) {
    TestUtils.deleteTestOutputDirectory(testInfo);
  }



  @Test
  void requiredPosition(TestInfo testInfo) throws Exception {
    Result rez = testMatchNamedAlleles(s_env, sf_definitionFile, new TestVcfBuilder(testInfo)
        .reference("NAT2")
        .missing("NAT2", "rs1801279")
        .generate());

    assertEquals(1, rez.getGeneCalls().size());
    GeneCall gc = rez.getGeneCalls().get(0);
    assertEquals("NAT2", gc.getGene());
    assertEquals(0, gc.getDiplotypes().size());
  }


  @Test
  void nat2Sheet(TestInfo testInfo) throws Exception {
    // this file based on https://docs.google.com/spreadsheets/d/1XjvmKHP7pp2t1qoqLKByD4bGjvf2Swpe13DqD18Ws-w/
    Path tsvFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/NAT2/tests.tsv");
    try (BufferedReader reader = Files.newBufferedReader(tsvFile)) {
      // header row
      String line = reader.readLine();
      assertNotNull(line);
      int row = 1;
      while ((line = reader.readLine()) != null) {
        row += 1;

        if (row == 33) {
          //continue;
        }

        String name = String.format("row_%03d", row);
        System.out.println("\n" + name);
        String[] data = line.split("\t");
        assertTrue(data.length > 7, name + " only has " + data.length + " columns, expecting at least 8");

        List<String> missingPositions = new ArrayList<>();
        if (!StringUtils.isBlank(data[0]) && !data[0].equals("none")) {
          String[] rsids = data[0].split("\\),");
          for (String rsid : rsids) {
            missingPositions.add(StringUtils.stripToNull(rsid.substring(0, rsid.indexOf("("))));
          }
        }

        boolean assumeReference = data[1].equals("reference");
        if (!assumeReference) {
          throw new UnsupportedOperationException("Currently only support assume-reference");
        }
        String[] rs1801279 = data[2].split("/");
        String[] rs1801280 = data[3].split("/");
        String[] rs1799930 = data[4].split("/");
        String[] rs1208 = data[5].split("/");
        String[] rs1799931 = data[6].split("/");

        TestVcfBuilder vcfBuilder = new TestVcfBuilder(testInfo, name)
            .variation("NAT2", "rs1801279", rs1801279[0], rs1801279[1])
            .variation("NAT2", "rs1801280", rs1801280[0], rs1801280[1])
            .variation("NAT2", "rs1799930", rs1799930[0], rs1799930[1])
            .variation("NAT2", "rs1208", rs1208[0], rs1208[1])
            .variation("NAT2", "rs1799931", rs1799931[0], rs1799931[1]);
        for (String missingPosition : missingPositions) {
          vcfBuilder.missing("NAT2", missingPosition);
        }

        Result rez = testMatchNamedAlleles(s_env, sf_definitionFile, vcfBuilder.generate(), true);
        GeneCall gc = rez.getGeneCalls().get(0);
        assertEquals("NAT2", gc.getGene());
        if (data[7].equals("not called")) {
          assertEquals(0, gc.getDiplotypes().size());
          continue;
        }
        List<String> matches = printMatches(gc);
        if (data.length > 8) {
          assertEquals(1, matches.size());
          assertEquals(List.of(data[8]), matches);
          continue;
        }

        List<String> output = sf_orSplitter.splitToList(data[7]);
        assertEquals(output.size(), matches.size(), name);
        assertEquals(output, matches);
      }
    }
  }
}
