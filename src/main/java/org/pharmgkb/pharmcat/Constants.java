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
  private static final List<String> VARIANT_GENES = List.of("ABCG2", "CACNA1S", "DPYD", "INFL3", "MT-RNR1", "RYR1",
      "VKORC1");
  public static final Set<String> PREFER_OUTSIDE_CALL = ImmutableSet.of("CYP2D6", "HLA-A", "HLA-B", "MT-RNR1");


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

  /**
   * Are the "alleles" specified for this gene variants (e.g. 1234A>G) instead of haplotypes (e.g. *1, *3)
   * <p>Note: this data should eventually come from PharmGKB and not be hard-coded</p>
   * @param geneSymbol a gene symbol
   * @return true if this genes alleles represent individual variants
   */
  public static boolean isVariantGene(String geneSymbol) {
    return VARIANT_GENES.contains(geneSymbol);
  }
}
