package org.pharmgkb.pharmcat.reporter.caller;

import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.phenotype.model.OutsideCall;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;
import org.pharmgkb.pharmcat.reporter.model.result.GeneReport;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;


/**
 * This is a JUnit test for {@link Cyp2d6CopyNumberCaller}.
 *
 * @author Mark Woon
 */
class Cyp2d6CopyNumberCallerTest {
  private static Env s_env;

  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  @Test
  void testWrongGene() {
    OutsideCall outsideCall = new OutsideCall(s_env, "ABCDE\t*1/*1x7", 1);

    GeneReport report = new GeneReport(outsideCall, s_env, DataSource.CPIC);
    assertEquals("ABCDE", report.getGene());

    Diplotype diplotype = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    assertEquals("ABCDE", outsideCall.getGene());

    assertThrows(IllegalArgumentException.class, () -> {
      Cyp2d6CopyNumberCaller.inferDiplotype(report, diplotype, s_env, DataSource.CPIC);
    });
  }


  @Test
  void testInferred() {
    OutsideCall outsideCall = new OutsideCall(s_env, "CYP2D6\t*1/*1x7", 1);

    GeneReport report = new GeneReport(outsideCall, s_env, DataSource.CPIC);
    assertEquals("CYP2D6", report.getGene());

    Diplotype diplotype = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    assertEquals("CYP2D6", outsideCall.getGene());

    Diplotype inferred = Cyp2d6CopyNumberCaller.inferDiplotype(report, diplotype, s_env, DataSource.CPIC);
    assertEquals("*1/*1x" + TextConstants.GTE + "3", inferred.getLabel());
  }




  @Test
  void testNoInfer1() {
    OutsideCall outsideCall = new OutsideCall(s_env, "CYP2D6\t*1/*42x7", 1);

    GeneReport report = new GeneReport(outsideCall, s_env, DataSource.CPIC);
    assertEquals("CYP2D6", report.getGene());

    Diplotype diplotype = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    assertEquals("CYP2D6", outsideCall.getGene());

    Diplotype inferred = Cyp2d6CopyNumberCaller.inferDiplotype(report, diplotype, s_env, DataSource.CPIC);
    assertEquals("*1/*42x7" , inferred.getLabel());
  }


  @Test
  void testNoInfer2() {
    OutsideCall outsideCall = new OutsideCall(s_env, "CYP2D6\t*1/*41x7", 1);

    GeneReport report = new GeneReport(outsideCall, s_env, DataSource.CPIC);
    assertEquals("CYP2D6", report.getGene());

    Diplotype diplotype = new Diplotype(outsideCall, s_env, DataSource.CPIC);
    assertEquals("CYP2D6", outsideCall.getGene());

    Diplotype inferred = Cyp2d6CopyNumberCaller.inferDiplotype(report, diplotype, s_env, DataSource.CPIC);
    assertEquals("*1/*41x7" , inferred.getLabel());
  }
}
