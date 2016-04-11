package org.cpic;

import org.junit.Test;
import org.pharmgkb.parser.vcf.VcfParser;

import java.io.BufferedReader;
import java.io.InputStreamReader;

/*
 * Simple test case showing how to use vcf-parser.
 *
 * @author Mark Woon
 */
public class VcfParserTest {

  @Test
  public void testReadingVcf() {

    try (BufferedReader reader = new BufferedReader(new InputStreamReader(VcfParserTest.class.getResourceAsStream("cyp2c9.vcf")))) {
      new VcfParser.Builder()
        .fromReader(reader)
        .parseWith((metadata, position, sampleData) -> {
          System.out.println(position.getChromosome());
          System.out.println(position.getPosition());
          System.out.println(position.getIds());
          System.out.println(position.getRef());
          System.out.println(position.getAltBases());
          System.out.println(position.getQuality());
        })
        .build().parse();
    } catch (Exception ex) {
      ex.printStackTrace();
    }
  }
}
