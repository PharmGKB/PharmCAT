package org.pharmgkb.pharmcat.phenotype.model;

import java.util.regex.Matcher;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.BadOutsideCallException;

import static org.junit.jupiter.api.Assertions.*;


/**
 * This is a JUnit test for {@link OutsideCall}.
 *
 * @author Mark Woon
 */
class OutsideCallTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  @Test
  void checkCyp2d6SuballelePattern() {
    Matcher m = OutsideCall.CYP2D6_SUBALLELE_PATTERN.matcher("*1.0001");
    assertTrue(m.matches());
    assertEquals("*1", m.group(1));

    m = OutsideCall.CYP2D6_SUBALLELE_PATTERN.matcher("*2");
    assertFalse(m.matches());


    OutsideCall oc = new OutsideCall(s_env, "CYP2D6\t*2/*2.001", 1);
    assertEquals("*2/*2", oc.getDiplotype());
    assertEquals(1, oc.getWarnings().size());

    oc = new OutsideCall(s_env, "CYP2D6\t*2.001/*3.024", 1);
    assertEquals("*2/*3", oc.getDiplotype());
    assertEquals(2, oc.getWarnings().size());
  }

  @Test
  void checkHlaSuballelePattern() {
    Matcher m = OutsideCall.HLA_SUBALLELE_PATTERN.matcher("*07:02:04:03");
    assertTrue(m.matches());
    assertEquals("*07:02", m.group(1));

    m = OutsideCall.HLA_SUBALLELE_PATTERN.matcher("*02:01");
    assertFalse(m.matches());

    OutsideCall oc = new OutsideCall(s_env, "HLA-A\t*02:02:02/*02:01", 1);
    assertEquals("*02:02/*02:01", oc.getDiplotype());
    assertEquals(1, oc.getWarnings().size());

    oc = new OutsideCall(s_env, "HLA-A\tA*02:02:02/A*02:01", 1);
    assertEquals("*02:02/*02:01", oc.getDiplotype());
    assertEquals(2, oc.getWarnings().size());

    oc = new OutsideCall(s_env, "HLA-A\tA*02:02:02/*02:01", 1);
    assertEquals("*02:02/*02:01", oc.getDiplotype());
    assertEquals(1, oc.getWarnings().size());

    assertThrows(BadOutsideCallException.class, () -> {
      new OutsideCall(s_env, "HLA-A\tB*02:02/*02:01", 1);
    }, "Invalid HLA-A allele: 'B*02:02'.");
  }


  @Test
  void testBad() {
    assertThrows(BadOutsideCallException.class, () -> {
      try {
        new OutsideCall(s_env, "CYP2D6", 1);
      } catch (BadOutsideCallException ex) {
        System.out.println(ex.getMessage());
        throw ex;
      }
    });

    assertThrows(BadOutsideCallException.class, () -> {
      try {
        new OutsideCall(s_env, "CYP2D6\t", 2);
      } catch (BadOutsideCallException ex) {
        System.out.println(ex.getMessage());
        throw ex;
      }
    });

    assertThrows(BadOutsideCallException.class, () -> {
      try {
        new OutsideCall(s_env, "CYP2D6\t   ", 3);
      } catch (BadOutsideCallException ex) {
        System.out.println(ex.getMessage());
        throw ex;
      }
    });

    assertThrows(BadOutsideCallException.class, () -> {
      try {
        new OutsideCall(s_env, "CYP2D6\tCYP2D6 *1/CYP2D6 *2/CYP2D6 *3", 4);
      } catch (BadOutsideCallException ex) {
        System.out.println(ex.getMessage());
        throw ex;
      }
    });
  }

  @Test
  void testGood() {
    new OutsideCall(s_env, "CYP2D6\t*1", 1);
    new OutsideCall(s_env, "CYP2D6\t*1\t",2);
    new OutsideCall(s_env, "CYP2D6\t*1\t\t",3);
    new OutsideCall(s_env, "CYP2D6\t\tIM", 4);
    new OutsideCall(s_env, "CYP2D6\t\t\t2.0", 5);
    new OutsideCall(s_env, "CYP2D6\t*1\tIM\t", 6);
    new OutsideCall(s_env, "CYP2D6\t*1\tIM\t2.0", 7);
    new OutsideCall(s_env, "CYP2D6\t*1\t\t2.0", 8);
    new OutsideCall(s_env, "CYP2D6\t\tIM\2.0", 9);

    OutsideCall oc = new OutsideCall(s_env, "UGT1A1\t*80+*28", 11);
    assertEquals(0, oc.getWarnings().size());

    oc = new OutsideCall(s_env, "UGT1A1\t*2/*80 + *28", 11);
    oc.getWarnings().forEach(System.out::println);
    assertEquals("*2/*80 + *28", oc.getDiplotype());
    assertEquals(2, oc.getWarnings().size());
    assertTrue(oc.getWarnings().stream().allMatch(w -> w.startsWith("Undocumented")));

    oc = new OutsideCall(s_env, "UGT1A1\t[*80 + *28]/*1", 11);
    //oc.getWarnings().forEach(System.out::println);
    assertEquals("[*80 + *28]/*1", oc.getDiplotype());
    assertEquals(1, oc.getWarnings().size());
    assertTrue(oc.getWarnings().get(0).startsWith("Undocumented"));
  }
}
