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
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.Haplotype;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.pharmgkb.pharmcat.reporter.model.DataSource.CPIC;
import static org.pharmgkb.pharmcat.reporter.model.DataSource.DPWG;


/**
 * This is a JUnit test for {@link DpydCaller}.
 *
 * @author Mark Woon
 */
class DpydCallerTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  @Test
  void infer_outsideCall_noMatch() {
    List<Diplotype> dips = DpydCaller.inferFromOutsideCall("*1/*5", s_env, CPIC);

    assertEquals(1, dips.size());
    Diplotype dip = dips.get(0);
    assertEquals("*1", dip.getAllele1().getName());
    assertEquals("*5", dip.getAllele2().getName());
    assertEquals(1, dip.getPhenotypes().size());
    Assertions.assertEquals(TextConstants.NA, dip.getPhenotypes().get(0));
  }


  @Test
  void inferPhased() {
    //  c.498G>A    - normal
    //  c.2582A>G   - normal
    String diplotype = "c.498G>A/c.2582A>G";
    checkInferred(diplotype, DataSource.CPIC, "c.498G>A", "c.2582A>G",
        "Normal function", "Normal function");
    checkInferred(diplotype, DPWG, "Reference", "Reference",
        "Normal function", "Normal function");

    //  c.2582A>G   - normal
    //  c.2846A>T   - decreased, DPWG
    diplotype = "c.2582A>G/c.2846A>T";
    checkInferred(diplotype, DataSource.CPIC, "c.2582A>G", "c.2846A>T",
        "Normal function", "Decreased function");
    checkInferred(diplotype, DPWG, "Reference", "c.2846A>T",
        "Normal function", "Decreased function");

    //  c.2582A>G   - normal
    //  c.2933A>G   - no function
    diplotype = "c.2582A>G/c.2933A>G";
    checkInferred(diplotype, DataSource.CPIC, "c.2582A>G", "c.2933A>G",
        "Normal function", "No function");
    checkInferred(diplotype, DPWG, "Reference", "c.2933A>G",
        "Normal function", GenePhenotype.UNASSIGNED_FUNCTION);

    //  c.2933A>G   - no function
    diplotype = "c.2933A>G/c.2933A>G";
    checkInferred(diplotype, DataSource.CPIC, "c.2933A>G", "c.2933A>G",
        "No function", "No function");
    checkInferred(diplotype, DPWG, "c.2933A>G", "c.2933A>G",
        GenePhenotype.UNASSIGNED_FUNCTION, GenePhenotype.UNASSIGNED_FUNCTION);

    //  c.498G>A    - normal
    //  c.2933A>G   - no function
    //  c.1905+1G>A (*2A) - no function, DPWG
    diplotype = "c.2933A>G + c.1905+1G>A (*2A)/c.498G>A";
    checkInferred(diplotype, DataSource.CPIC, "c.1905+1G>A (*2A)", "c.498G>A",
        "No function", "Normal function");
    checkInferred(diplotype, DPWG, "c.1905+1G>A (*2A)", "Reference",
        "No function", "Normal function");

    //  c.498G>A    - normal
    //  c.2933A>G   - no function
    //  c.1905+1G>A (*2A) - no function, DPWG
    diplotype = "c.498G>A + c.2933A>G/c.1905+1G>A (*2A)";
    checkInferred(diplotype, DataSource.CPIC, "c.2933A>G", "c.1905+1G>A (*2A)",
        "No function", "No function");
    checkInferred(diplotype, DPWG, "c.2933A>G", "c.1905+1G>A (*2A)",
        GenePhenotype.UNASSIGNED_FUNCTION, "No function");

    //  c.498G>A    - normal
    //  c.2582A>G   - normal
    //  c.2846A>T   - decreased, DPWG
    //  c.2933A>G   - no function
    diplotype = "c.498G>A + c.2582A>G/c.2846A>T + c.2933A>G";
    checkInferred(diplotype, DataSource.CPIC, "c.498G>A", "c.2933A>G",
        "Normal function", "No function");
    checkInferred(diplotype, DPWG, "Reference", "c.2933A>G",
        "Normal function", GenePhenotype.UNASSIGNED_FUNCTION);
  }

  private void checkInferred(String diplotype, DataSource source, String a1, String a2, String f1, String f2) {
    List<Diplotype> dips = DpydCaller.inferFromOutsideCall(diplotype, s_env, source);

    assertEquals(1, dips.size());
    Diplotype dip = dips.get(0);
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

    List<Diplotype> dips = DpydCaller.inferFromHaplotypeMatches(matches, s_env, CPIC);
    Diplotype dip = dips.get(0);
    assertEquals(na4.getName(), dip.getAllele1().getName());
    assertEquals(na3.getName(), dip.getAllele2().getName());
    assertEquals("No function", dip.getAllele1().getFunction());
    assertEquals("Decreased function", dip.getAllele2().getFunction());

    dips = DpydCaller.inferFromHaplotypeMatches(matches, s_env, DPWG);
    dip = dips.get(0);
    assertEquals(na4.getName(), dip.getAllele1().getName());
    assertEquals(na3.getName(), dip.getAllele2().getName());
    assertEquals(GenePhenotype.UNASSIGNED_FUNCTION, dip.getAllele1().getFunction());
    assertEquals("Decreased function", dip.getAllele2().getFunction());
  }


  @Test
  void testComparator() {
    DpydCaller.DpydActivityComparator comparator = new DpydCaller.DpydActivityComparator(s_env);

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
}
