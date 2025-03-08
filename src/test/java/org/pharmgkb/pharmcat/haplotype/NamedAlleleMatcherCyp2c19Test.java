package org.pharmgkb.pharmcat.haplotype;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import com.google.common.collect.Lists;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestInfo;
import org.pharmgkb.common.util.PathUtils;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.TestUtils;
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;

import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.assertDiplotypePairs;
import static org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcherTest.testMatchNamedAlleles;


/**
 * CYP2C19-specific tests for {@link NamedAlleleMatcher}.
 *
 * @author Mark Woon
 */
public class NamedAlleleMatcherCyp2c19Test {
  private static final Path sf_definitionFile = DataManager.DEFAULT_DEFINITION_DIR.resolve("CYP2C19_translation.json");
  private static Env s_env = null;

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
  void s1s1(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("*1/*1", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*1")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .generate()));
  }


  @Test
  void s1s2(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("*1/*2", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*1")
            .variation("CYP2C19", "rs12769205", "A", "G")
            .variation("CYP2C19", "rs4244285", "G", "A")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .generate()));
  }


  @Test
  void s1s4(TestInfo testInfo) throws Exception {
    // Test top condidate is *1/*4 (over *4/*17)
    assertDiplotypePairs("*1/*4", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*4")
            .variation("CYP2C19", "rs12248560", "C", "T")
            .variation("CYP2C19", "rs28399504", "A", "G")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .generate(),
        true));

  }

  @Test
  void s1s17s1s4bMissingCalls(TestInfo testInfo) throws Exception {
    // Test *1/*4 - but with 94762706	rs28399504	A	G	.	PASS	star-4a-4b	GT	./. to test partial call
    assertDiplotypePairs(Lists.newArrayList("*1/*4", "*1/*17"), testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*4 and *1/*17")
            .variation("CYP2C19", "rs12248560", "C", "T")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .missing("CYP2C19", "rs28399504")
            .generate(),
        true));
  }


  @Test
  void s1s17s1s4bMissingMoreCalls(TestInfo testInfo) throws Exception {
    // Test *1/*4 - but many calls deleted from the vcf.
    // TODO: returned results are correct.  In the html report *4 A and *15 are reported as excluded.
    //  However as far as I can tell all positions for *28 are also excluded, so why doesn't this make the list?
    //  Presume this is becuase it's multi position, so similar to *1, but may be worth discussing.
    assertDiplotypePairs(Lists.newArrayList("*1/*4", "*1/*17"), testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*4 and *1/*17")
            .variation("CYP2C19", "rs12248560", "C", "T")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .missing("CYP2C19", "rs28399504")
            .missing("CYP2C19", "rs17882687")
            .missing("CYP2C19", "rs1564656981")
            .missing("CYP2C19", "rs1564657013")
            .missing("CYP2C19", "rs1564660997")
            .missing("CYP2C19", "rs1288601658")
            .missing("CYP2C19", "rs17885179")
            .missing("CYP2C19", "rs375781227")
            .missing("CYP2C19", "rs72558186")
            .missing("CYP2C19", "rs113934938")
            .generate(),
        true));
  }


  @Test
  void s1s28(TestInfo testInfo) throws Exception {
    // Test *1 *28. The longest possible match should win.
    assertDiplotypePairs("*1/*28", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*28")
            .variation("CYP2C19", "rs17882687", "A", "C")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .variation("CYP2C19", "rs113934938", "G", "A")
            .generate()));
  }


  @Test
  void s2s2(TestInfo testInfo) throws Exception {
    // Test simple case of one homozygous snp
    // TODO: description of test no longer matches data (have to add missing position for this to pass)
    //  should we update data to specify missing position(s) instead?
    assertDiplotypePairs("*2/*2", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*2/*2")
            .variation("CYP2C19", "rs12769205", "G", "G")
            .variation("CYP2C19", "rs4244285", "A", "A")
            .missing("CYP2C19", "rs3758581")
            .generate()));
  }


  @Test
  void s2s3(TestInfo testInfo) throws Exception {
    // TODO: description of test no longer matches data (have to add missing position for this to pass)
    //  should we update data to specify missing position(s) instead?
    assertDiplotypePairs("*2/*3", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*2/*3")
            .variation("CYP2C19", "rs12769205", "A", "G")
            .variation("CYP2C19", "rs4986893", "G", "A")
            .variation("CYP2C19", "rs4244285", "G", "A")
            .missing("CYP2C19", "rs3758581")
            .generate()));
  }


  @Test
  void s4as4b(TestInfo testInfo) throws Exception {
    // Test *4a *4b. het and homo rsids
    // TODO: description of test no longer matches data (have to add missing position for this to pass)
    //  should we update data to specify missing position(s) instead?
    assertDiplotypePairs("*4/*4", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*4/*4")
            .variation("CYP2C19", "rs12248560", "C", "T")
            .variation("CYP2C19", "rs28399504", "G", "G")
            .missing("CYP2C19", "rs3758581")
            .generate()));
  }


  @Test
  void s4bs17(TestInfo testInfo) throws Exception {
    // Test *4b/*17
    // TODO: description of test no longer matches data (have to add missing position for this to pass)
    //  should we update data to specify missing position(s) instead?
    assertDiplotypePairs("*4/*17", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*4/*17")
            .variation("CYP2C19", "rs12248560", "T", "T")
            .variation("CYP2C19", "rs28399504", "A", "G")
            .missing("CYP2C19", "rs3758581")
            .generate()));
  }


  @Test
  void s15s28(TestInfo testInfo) throws Exception {
    // Test *15 *28. The shared position is homo
    // TODO: description of test no longer matches data (have to add missing position for this to pass)
    //  should we update data to specify missing position(s) instead?
    assertDiplotypePairs(Lists.newArrayList("*15/*28", "*28/*39"), testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*15/*28 and *28/*39")
            .variation("CYP2C19", "rs17882687", "C", "C")
            .variation("CYP2C19", "rs113934938", "G", "A")
            .missing("CYP2C19", "rs3758581")
            .missing("CYP2C19", "rs17885179")
            .generate()));
  }


  @Test
  void sUnks17(TestInfo testInfo) throws Exception {
    // Test *Unk/*17 - only one haplotype matches, so no diploid match
    assertDiplotypePairs(Lists.newArrayList(), testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "unknown/*17")
            .variation("CYP2C19", "rs12248560", "T", "T")
            .variation("CYP2C19", "rs28399504", "A", "G")
            .generate()));
  }


  @Test
  void rs12769205missingrs4244285het(TestInfo testInfo) throws Exception {
    assertDiplotypePairs(Lists.newArrayList("*1/*2", "*2/*35"), testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*2 and *2/*35")
            .variation("CYP2C19", "rs4244285", "G", "A")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .missing("CYP2C19", "rs12769205")
            .generate()));
  }


  // removing rs12769205hetrs4244285missing test since rs12769205call duplicates it with new 2C19 definitions

  @Test
  void rs12769205call(TestInfo testInfo) throws Exception {
    // Test no call, but reporter can give output based on rs12769205
    // TODO: description of test no longer matches data
    assertDiplotypePairs("*1/*35", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*35")
            .variation("CYP2C19", "rs12769205", "A", "G")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .generate()));
  }

  @Test
  void s4bs17rs28399504missing(TestInfo testInfo) throws Exception {
    // rs28399504 missing
    assertDiplotypePairs(Lists.newArrayList(), testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "no call")
            .variation("CYP2C19", "rs12248560", "T", "T")
            .missing("CYP2C19", "rs28399504")
            .generate()));

    Path vcfFile = PathUtils.getPathToResource("org/pharmgkb/pharmcat/haplotype/cyp2c19/s4bs17rs28399504missing.vcf");
    List<String> expectedMatches = Lists.newArrayList("*4/*4", "*4/*17", "*17/*17");

    Result result = testMatchNamedAlleles(s_env, sf_definitionFile, vcfFile);
    assertDiplotypePairs(expectedMatches, result);
  }


  @Test
  void s1s1rs12248560missing(TestInfo testInfo) throws Exception {
    // rs12248560 missing
    assertDiplotypePairs(Lists.newArrayList("*1/*1", "*1/*17", "*17/*17"), testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*1/*1, *1/*17, *17/*17")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .missing("CYP2C19", "rs12248560")
            .generate()));
  }


  @Test
  void s4s17het(TestInfo testInfo) throws Exception {
    assertDiplotypePairs("*4/*17", testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "*4/*17")
            .variation("CYP2C19", "rs12248560", "T", "T")
            .variation("CYP2C19", "rs28399504", "A", "G")
            .variation("CYP2C19", "rs3758581", "G", "G")
            .generate()));
  }


  @Test
  void s2s35(TestInfo testInfo) throws Exception {
    assertDiplotypePairs(Lists.newArrayList(), testMatchNamedAlleles(s_env, sf_definitionFile,
        new TestVcfBuilder(testInfo, "no call")
            .variation("CYP2C19", "rs12769205", "A", "G")
            .variation("CYP2C19", "rs4244285", "A", "A")
            .generate()));
  }
}
