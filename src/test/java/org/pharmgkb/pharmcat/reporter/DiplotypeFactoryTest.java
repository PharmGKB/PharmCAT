package org.pharmgkb.pharmcat.reporter;

import java.util.ArrayList;
import java.util.List;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.reporter.model.DataSource;
import org.pharmgkb.pharmcat.reporter.model.result.Diplotype;

import static org.junit.jupiter.api.Assertions.assertEquals;


/**
 * This is a JUnit test for {@link DiplotypeFactory}.
 *
 * @author Mark Woon
 */
class DiplotypeFactoryTest {
  private static Env s_env;


  @BeforeAll
  static void prepare() throws Exception {
    s_env = new Env();
  }


  @Test
  void makeLeastFunctionDiplotypes_outsideCall() {
    DiplotypeFactory df = new DiplotypeFactory("DPYD", s_env);
    List<String> matches = new ArrayList<>();
    matches.add("*1/*5");
    List<Diplotype> dips = df.makeLeastFunctionDiplotypes(matches, DataSource.CPIC, true);

    assertEquals(1, dips.size());
    Diplotype dip = dips.get(0);
    assertEquals("*1", dip.getAllele1().getName());
    assertEquals("*5", dip.getAllele2().getName());
  }
}
