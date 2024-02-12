package org.pharmgkb.pharmcat;

import java.util.List;
import java.util.Set;
import com.google.common.collect.ImmutableSet;
import org.apache.commons.lang3.StringUtils;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;


/**
 * This class contains constants for use in PharmCAT.
 *
 * @author Mark Woon
 */
public class Constants {
  private static final List<String> LOWEST_FUNCTION_GENES = List.of("DPYD", "RYR1");
  /**
   * Genes that are on the X, Y, or M chromosome.
   * <p>
   * These genes have custom callers that can potentially infer a diplotype even though the {@link NamedAlleleMatcher}
   * cannot make a call.
   */
  private static final Set<String> SINGLE_PLOIDY = ImmutableSet.of("G6PD", "MT-RNR1");


  /**
   * Does this gene use the lowest-function algorithm for determining recommendations?
   */
  public static boolean isLowestFunctionGene(String gene) {
    return LOWEST_FUNCTION_GENES.contains(gene);
  }


  /**
   * Is the specified gene single ploidy (i.e. it is on the X, Y, or M chromosome that might not occur in a pair)?.
   *
   * @param gene a gene symbol
   * @return true if the gene occurs on a single chromosome, not a pair
   */
  public static boolean isSinglePloidy(String gene) {
    return StringUtils.isNotBlank(gene) && SINGLE_PLOIDY.contains(gene);
  }
}
