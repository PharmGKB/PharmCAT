package org.pharmgkb.pharmcat;

import java.util.List;


/**
 * This class contains constants for use in PharmCAT.
 *
 * @author Mark Woon
 */
public class Constants {
  private static final List<String> LOWEST_FUNCTION_GENES = List.of("DPYD", "RYR1");


  public static boolean isLowestFunctionGene(String gene) {
    return LOWEST_FUNCTION_GENES.contains(gene);
  }
}
