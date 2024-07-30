package org.pharmgkb.pharmcat.reporter.caller;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.phenotype.model.GenePhenotype;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.pharmgkb.pharmcat.reporter.model.DataSource.CPIC;
import static org.pharmgkb.pharmcat.reporter.model.DataSource.DPWG;


/**
 * This is a JUnit test for {@link LowestFunctionGeneCaller}.
 *
 * @author Mark Woon
 */
class LowestFunctionGeneCallerTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  @Test
  void infer_outsideCall_noMatch() {
    OutsideCall call = new OutsideCall(s_env, "DPYD\t*1/*5", 1);
    List<Diplotype> dips = LowestFunctionGeneCaller.inferFromOutsideCall(call, s_env, CPIC);

    assertEquals(1, dips.size());
    Diplotype dip = dips.get(0);
    assertNotNull(dip.getAllele1());
    assertEquals("*1", dip.getAllele1().getName());
    assertNotNull(dip.getAllele2());
    assertEquals("*5", dip.getAllele2().getName());
    assertEquals(1, dip.getPhenotypes().size());
    Assertions.assertEquals(TextConstants.NA, dip.getPhenotypes().get(0));
  }


  @Test
  void inferPhased() {
    //  Reference   - normal
    //  c.2846A>T   - decreased, DPWG
    String diplotype = "Reference/c.2846A>T";
    checkInferred(diplotype, DataSource.CPIC, "Reference", "c.2846A>T",
        "Normal function", "Decreased function");
    checkInferred(diplotype, DPWG, "Reference", "c.2846A>T",
        "Normal function", "Decreased function");

    //  Reference   - normal
    //  c.2933A>G   - no function
    diplotype = "Reference/c.2933A>G";
    checkInferred(diplotype, DataSource.CPIC, "Reference", "c.2933A>G",
        "Normal function", "No function");
    checkInferred(diplotype, DPWG, "Reference", "c.2933A>G",
        "Normal function", GenePhenotype.UNASSIGNED_FUNCTION);

    //  c.2933A>G   - no function
    diplotype = "c.2933A>G/c.2933A>G";
    checkInferred(diplotype, DataSource.CPIC, "c.2933A>G", "c.2933A>G",
        "No function", "No function");
    checkInferred(diplotype, DPWG, "c.2933A>G", "c.2933A>G",
        GenePhenotype.UNASSIGNED_FUNCTION, GenePhenotype.UNASSIGNED_FUNCTION);

    //  Reference    - normal
    //  c.2933A>G   - no function
    //  c.1905+1G>A (*2A) - no function, DPWG
    diplotype = "Reference/[c.2933A>G + c.1905+1G>A (*2A)]";
    checkInferred(diplotype, DataSource.CPIC, "Reference","c.1905+1G>A (*2A)",
        "Normal function", "No function");
    checkInferred(diplotype, DPWG, "Reference", "c.1905+1G>A (*2A)",
        "Normal function", "No function");

    //  Reference    - normal
    //  c.2933A>G   - no function
    //  c.1905+1G>A (*2A) - no function, DPWG
    diplotype = "c.1905+1G>A (*2A)/[Reference + c.2933A>G]";
    checkInferred(diplotype, DataSource.CPIC, "c.1905+1G>A (*2A)", "c.2933A>G",
        "No function", "No function");
    checkInferred(diplotype, DPWG, "c.1905+1G>A (*2A)", "c.2933A>G",
        "No function", GenePhenotype.UNASSIGNED_FUNCTION);

    //  Reference    - normal
    //  c.2846A>T   - decreased, DPWG
    //  c.2933A>G   - no function
    diplotype = "Reference/[c.2846A>T + c.2933A>G]";
    checkInferred(diplotype, DataSource.CPIC, "Reference", "c.2933A>G",
        "Normal function", "No function");
    checkInferred(diplotype, DPWG, "Reference", "c.2933A>G",
        "Normal function", GenePhenotype.UNASSIGNED_FUNCTION);
  }

  private void checkInferred(String diplotype, DataSource source, String a1, String a2, String f1, String f2) {
    OutsideCall call = new OutsideCall(s_env, "DPYD\t" + diplotype, 1);
    List<Diplotype> dips = LowestFunctionGeneCaller.inferFromOutsideCall(call, s_env, source);

    assertEquals(1, dips.size());
    Diplotype dip = dips.get(0);
    assertNotNull(dip.getAllele1());
    assertNotNull(dip.getAllele2());
    assertEquals(a1, dip.getAllele1().getName());
    assertEquals(f1, dip.getAllele1().getFunction());
    assertEquals(a2, dip.getAllele2().getName());
    assertEquals(f2, dip.getAllele2().getFunction());
  }


  @Test
  void inferDiplotypes_matcherUnphased() {
    //  Reference    - normal
    //  c.2846A>T   - decreased, DPWG
    //  c.2933A>G   - no function
    NamedAllele na1 = new NamedAllele("1", "Reference", new String[0], new String[0], false);
    NamedAllele na3 = new NamedAllele("3", "c.2846A>T", new String[0], new String[0], false);
    NamedAllele na4 = new NamedAllele("4", "c.2933A>G", new String[0], new String[0], false);

    List<HaplotypeMatch> matches = new ArrayList<>();
    matches.add(new HaplotypeMatch(na1));
    matches.add(new HaplotypeMatch(na3));
    matches.add(new HaplotypeMatch(na4));

    DiplotypeFactory diplotypeFactory = new DiplotypeFactory("DPYD", s_env);

    List<Diplotype> dips = LowestFunctionGeneCaller.inferFromHaplotypeMatches("DPYD", s_env, CPIC,
        diplotypeFactory, matches);
    Diplotype dip = dips.get(0);
    assertNotNull(dip.getAllele1());
    assertNotNull(dip.getAllele2());
    assertEquals(na4.getName(), dip.getAllele1().getName());
    assertEquals(na3.getName(), dip.getAllele2().getName());
    assertEquals("No function", dip.getAllele1().getFunction());
    assertEquals("Decreased function", dip.getAllele2().getFunction());

    dips = LowestFunctionGeneCaller.inferFromHaplotypeMatches("DPYD", s_env, DPWG, diplotypeFactory, matches);
    dip = dips.get(0);
    assertNotNull(dip.getAllele1());
    assertNotNull(dip.getAllele2());
    assertEquals(na4.getName(), dip.getAllele1().getName());
    assertEquals(na3.getName(), dip.getAllele2().getName());
    assertEquals(GenePhenotype.UNASSIGNED_FUNCTION, dip.getAllele1().getFunction());
    assertEquals("Decreased function", dip.getAllele2().getFunction());
  }


  @Test
  void testComparator() {
    LowestFunctionGeneCaller.DpydActivityComparator comparator = new LowestFunctionGeneCaller.DpydActivityComparator(s_env);

    Haplotype hapRef = new Haplotype("DPYD", "Reference");  // normal, DPWG
    Haplotype hap2846 = new Haplotype("DPYD", "c.2846A>T");  // decreased, DPWG
    Haplotype hap2933 = new Haplotype("DPYD", "c.2933A>G");  // no function

    // decreased before normal
    assertEquals(-1, comparator.compare(hap2846, hapRef));
    // no func before decreased
    assertEquals(1, comparator.compare(hap2846, hap2933));
    // no func before normal
    assertEquals(-1, comparator.compare(hap2933, hapRef));
  }
}
