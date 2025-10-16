package org.pharmgkb.pharmcat.reporter.caller;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.DiplotypeFactory;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;


/**
 * This is a JUnit test for {@link LowestFunctionGeneCaller}.
 *
 * @author Mark Woon
 */
class LowestFunctionGeneCallerTest {
  private static Env s_env;
  private static DiplotypeFactory s_diplotypeFactory;


  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
    s_diplotypeFactory = new DiplotypeFactory("DPYD", s_env);
  }


  @Test
  void infer_outsideCall_noMatch() {
    OutsideCall call = new OutsideCall(s_env, "DPYD\t*1/*5", 1);
    List<Diplotype> dips = LowestFunctionGeneCaller.inferFromOutsideCall(call, s_env);

    assertEquals(1, dips.size());
    Diplotype dip = dips.get(0);
    assertNotNull(dip.getAllele1());
    assertEquals("*1", dip.getAllele1().getName());
    assertNotNull(dip.getAllele2());
    assertEquals("*5", dip.getAllele2().getName());
    assertEquals(1, dip.getPhenotypes().size());
    Assertions.assertEquals(TextConstants.INDETERMINATE, dip.getPhenotypes().get(0));
  }


  @Test
  void inferPhased() {
    //  c.498G>A    - normal
    //  c.2582A>G   - normal
    String diplotype = "c.498G>A/c.2582A>G";
    checkInferred(diplotype, "c.498G>A", "c.2582A>G",
        "Normal function", "Normal function");

    //  c.2582A>G   - normal
    //  c.2846A>T   - decreased, DPWG
    diplotype = "c.2582A>G/c.2846A>T";
    checkInferred(diplotype, "c.2582A>G", "c.2846A>T",
        "Normal function", "Decreased function");

    //  c.2582A>G   - normal
    //  c.2933A>G   - no function
    diplotype = "c.2582A>G/c.2933A>G";
    checkInferred(diplotype, "c.2582A>G", "c.2933A>G",
        "Normal function", "No function");

    //  c.2933A>G   - no function
    diplotype = "c.2933A>G/c.2933A>G";
    checkInferred(diplotype, "c.2933A>G", "c.2933A>G",
        "No function", "No function");

    //  c.498G>A    - normal
    //  c.2933A>G   - no function
    //  c.1905+1G>A (*2A) - no function, DPWG
    diplotype = "c.498G>A/[c.2933A>G + c.1905+1G>A (*2A)]";
    checkInferred(diplotype, "c.498G>A","c.1905+1G>A (*2A)",
        "Normal function", "No function");

    //  c.498G>A    - normal
    //  c.2933A>G   - no function
    //  c.1905+1G>A (*2A) - no function, DPWG
    diplotype = "c.1905+1G>A (*2A)/[c.498G>A + c.2933A>G]";
    checkInferred(diplotype, "c.1905+1G>A (*2A)", "c.2933A>G",
        "No function", "No function");

    //  c.498G>A    - normal
    //  c.2582A>G   - normal
    //  c.2846A>T   - decreased, DPWG
    //  c.2933A>G   - no function
    diplotype = "[c.498G>A + c.2582A>G]/[c.2846A>T + c.2933A>G]";
    checkInferred(diplotype, "c.498G>A", "c.2933A>G",
        "Normal function", "No function");
  }

  private void checkInferred(String diplotype, String a1, String a2, String f1, String f2) {
    OutsideCall call = new OutsideCall(s_env, "DPYD\t" + diplotype, 1);
    List<Diplotype> dips = LowestFunctionGeneCaller.inferFromOutsideCall(call, s_env);

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
    //  c.498G>A    - normal
    //  c.2582A>G   - normal
    //  c.2846A>T   - decreased, DPWG
    //  c.2933A>G   - no function
    NamedAllele na1 = new NamedAllele("1", "c.498G>A", new String[0], new String[0], false);
    NamedAllele na2 = new NamedAllele("2", "c.2582A>G", new String[0], new String[0], false);
    NamedAllele na3 = new NamedAllele("3", "c.2846A>T", new String[0], new String[0], false);
    NamedAllele na4 = new NamedAllele("4", "c.2933A>G", new String[0], new String[0], false);

    List<HaplotypeMatch> matches = new ArrayList<>();
    matches.add(new HaplotypeMatch(na1));
    matches.add(new HaplotypeMatch(na2));
    matches.add(new HaplotypeMatch(na3));
    matches.add(new HaplotypeMatch(na4));

    DiplotypeFactory diplotypeFactory = new DiplotypeFactory("DPYD", s_env);

    List<Diplotype> dips = LowestFunctionGeneCaller.inferFromHaplotypeMatches("DPYD", s_env,
        diplotypeFactory, matches);
    Diplotype dip = dips.get(0);
    assertNotNull(dip.getAllele1());
    assertNotNull(dip.getAllele2());
    assertEquals(na4.getName(), dip.getAllele1().getName());
    assertEquals(na3.getName(), dip.getAllele2().getName());
    assertEquals("No function", dip.getAllele1().getFunction());
    assertEquals("Decreased function", dip.getAllele2().getFunction());
  }


  private void checkInferredDpyd(List<DiplotypeMatch> matches, String cpicDip, String cpicActivityScore, String dpwgDip) {
    List<Diplotype> rez = LowestFunctionGeneCaller.inferFromDiplotypes("DPYD", s_env, s_diplotypeFactory,
        matches);
    assertEquals(1, rez.size());
    Diplotype finalDip = rez.get(0);
    assertEquals(cpicDip, finalDip.buildLabel(true));
    assertEquals(cpicActivityScore, finalDip.getActivityScore(), "Incorrect inferred CPIC diplotype");
  }


  @Test
  void inferDiplotypes_matcherPhaseSets() {

    // rs1801267 C > T - c.2657G>A (*9B) - Activity Value 1.0 (Normal function), not in DPWG
    NamedAllele na9b = new NamedAllele("1", "c.2657G>A (*9B)", new String[0], new String[0], false);
    HaplotypeMatch hm9b = new HaplotypeMatch(na9b);
    // rs59086055 G > A - c.1774C>T - Activity Value 0.0 (No function), not in DPWG
    NamedAllele na1774 = new NamedAllele("2", "c.1774C>T", new String[0], new String[0], false);
    HaplotypeMatch hm1774 = new HaplotypeMatch(na1774);
    // rs186169810 A > A - c.1314T>G - Activity Value 0.5 (Decreased function), not in DPWG
    NamedAllele na1314 = new NamedAllele("3", "c.1314T>G", new String[0], new String[0], false);

    // rs55886062 A > C - c.1679T>G (*13) - Activity Value 0.0 (No function)
    NamedAllele na13 = new NamedAllele("4", "c.1679T>G (*13)", new String[0], new String[0], false);
    HaplotypeMatch hm13 = new HaplotypeMatch(na13);
    // rs67376798 T > T - c.2846A>T - Activity Value 0.5 (Decreased function)
    NamedAllele na2846 = new NamedAllele("4", "c.2846A>T", new String[0], new String[0], false);


    // *9B/*13 -> normal (no DPWG)/no function
    List<DiplotypeMatch> matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm9b, hm13, null));

    checkInferredDpyd(matches, "c.1679T>G (*13)/c.2657G>A (*9B)", "1.0",
        "Reference/c.1679T>G (*13)");

    //-----
    // c.1774C>T/*13 -> no function (no DPWG)/no function
    matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm1774, hm13, null));

    checkInferredDpyd(matches, "c.1679T>G (*13)/c.1774C>T", "0.0",
        "c.1679T>G (*13)/c.1774C>T");


    //-----
    // *9B/[c.1774C>T + *13] -> normal (no DPWG)/no function (no DPWG) + no function
    // should pick DPWG no function over CPIC no function
    matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm9b, new CombinationMatch(new VariantLocus[0], "", List.of(na1774, na13), null), null));

    checkInferredDpyd(matches, "c.1679T>G (*13)/c.2657G>A (*9B)", "1.0",
        "Reference/c.1679T>G (*13)");


    //-----
    // *9B/[c.1774C>T + *13] -> normal (no DPWG)/no function (no DPWG) + no function
    // should pick DPWG no function over CPIC no function
    matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm9b, new CombinationMatch(new VariantLocus[0], "", List.of(na1774, na13), null), null));

    checkInferredDpyd(matches, "c.1679T>G (*13)/c.2657G>A (*9B)", "1.0",
        "Reference/c.1679T>G (*13)");


    //-----
    // *9B/[c.1774C>T + *13] -> normal (no DPWG)/no function (no DPWG) + no function
    // [*9B + c.1774C>T]/*13 -> no function (no DPWG) + normal (no DPWG)/no function
    // should pick DPWG no function over CPIC no function
    matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm9b, new CombinationMatch(new VariantLocus[0], "", List.of(na1774, na13), null), null));
    matches.add(new DiplotypeMatch(hm13, new CombinationMatch(new VariantLocus[0], "", List.of(na1774, na9b), null), null));

    checkInferredDpyd(matches, "c.1679T>G (*13)/c.1774C>T", "0.0",
        "c.1679T>G (*13)/c.1774C>T");
  }


  @Test
  void testDpydActivityComparator() {
    LowestFunctionGeneCaller.DpydActivityComparator comparator = new LowestFunctionGeneCaller.DpydActivityComparator(s_env);

    Haplotype hapRef = new Haplotype("DPYD", "Reference");  // normal, DPWG
    Haplotype hap498 = new Haplotype("DPYD", "c.498G>A");   // normal
    Haplotype hap2582 = new Haplotype("DPYD", "c.2582A>G"); // normal
    Haplotype hap2846 = new Haplotype("DPYD", "c.2846A>T");  // decreased, DPWG
    Haplotype hap2933 = new Haplotype("DPYD", "c.2933A>G");  // no function

    // decreased before normal
    assertEquals(-1, comparator.compare(hap2846, hapRef));
    // no func before decreased
    assertEquals(1, comparator.compare(hap2846, hap2933));
    // no func before normal
    assertEquals(-1, comparator.compare(hap2933, hapRef));

    // prefer DPWG
    assertEquals(1, comparator.compare(hap498, hapRef));
    assertEquals(-1, comparator.compare(hapRef, hap498));

    // order by name
    assertEquals(1, comparator.compare(hap2582, hap498));
    assertEquals(-1, comparator.compare(hap498, hap2582));
  }


  @Test
  void testRyr1DiplotypeComparator() {

    // c.38T>G - Malignant Hyperthermia associated
    NamedAllele na38 = new NamedAllele("1", "c.38T>G", new String[0], new String[0], false);
    HaplotypeMatch hm38 = new HaplotypeMatch(na38);
    // c.51_53del - Uncertain function
    NamedAllele na51 = new NamedAllele("2", "c.51_53del", new String[0], new String[0], false);
    HaplotypeMatch hm51 = new HaplotypeMatch(na51);


    // c.38T>G/c.51_53del - malignant/uncertain
    List<DiplotypeMatch> matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm38, hm51, null));
    checkInferredRyr1(matches, "c.38T>G/c.51_53del");

    // c.38T>G/c.51_53del - malignant/uncertain
    // c.51_53del/c.51_53del - uncertain/uncertain
    matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm38, hm51, null));
    matches.add(new DiplotypeMatch(hm51, hm51, null));
    checkInferredRyr1(matches, "c.38T>G/c.51_53del");

    // c.38T>G/c.38T>G - malignant/malignant
    // c.51_53del/c.51_53del - uncertain/uncertain
    matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm38, hm38, null));
    matches.add(new DiplotypeMatch(hm51, hm51, null));
    checkInferredRyr1(matches, "c.38T>G/c.38T>G");

    // c.38T>G/c.38T>G - malignant/malignant
    // c.38T>G/c.51_53del - malignant/uncertain
    matches = new ArrayList<>();
    matches.add(new DiplotypeMatch(hm38, hm38, null));
    matches.add(new DiplotypeMatch(hm38, hm51, null));
    checkInferredRyr1(matches, "c.38T>G/c.38T>G");
  }

  private void checkInferredRyr1(List<DiplotypeMatch> matches, String cpicDip) {
    List<Diplotype> rez = LowestFunctionGeneCaller.inferFromDiplotypes("RYR1", s_env, s_diplotypeFactory,
        matches);
    assertEquals(1, rez.size());
    Diplotype finalDip = rez.get(0);
    assertEquals(cpicDip, finalDip.buildLabel(true));
  }
}
