package org.pharmgkb.pharmcat;

import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.ImmutableSortedSet;
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
  private static final Set<String> X_CHROMOSOME_GENES = ImmutableSet.of("G6PD");

  /**
   * Genes that should be labeled as "Named Variants".
   */
  private static final Set<String> VARIANT_GENES = ImmutableSet.of(
      "ABCG2", "CACNA1S", "CFTR", "DPYD", "G6PD", "INFL3", "MT-RNR1", "RYR1", "VKORC1");

  public static final Set<String> ALLELE_PRESENCE_GENES = ImmutableSortedSet.of("HLA-A", "HLA-B");
  private static final Set<String> ACTIVITY_SCORE_GENES = ImmutableSortedSet.of("CYP2C9", "CYP2D6", "DPYD");

  public static final SortedSet<String> PREFER_OUTSIDE_CALL = ImmutableSortedSet.of("CYP2D6", "HLA-A", "HLA-B", "MT-RNR1");

  public static final SortedSet<String> HLA_A_ALLELES = ImmutableSortedSet.of("*31:01");
  public static final SortedSet<String> HLA_B_ALLELES = ImmutableSortedSet.of("*15:02", "*57:01", "*58:01");


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
   * Are the "alleles" specified for this gene variants (e.g. 1234A>G) instead of haplotypes (e.g. *1, *3)?
   * <p>
   * Note: this data should eventually come from PharmGKB and not be hard-coded
   *
   * @param geneSymbol a gene symbol
   * @return true if this genes alleles represent individual variants
   */
  public static boolean isVariantGene(String geneSymbol) {
    return VARIANT_GENES.contains(geneSymbol);
  }

  /**
   * Returns true if the gene does not use allele function to assign phenotype but instead relies on the presence or
   * absense of alleles for its phenotypes (e.g. HLA's).
   *
   * @param gene the gene symbol
   * @return true if this gene assigns phenotype based on allele presence
   */
  public static boolean isAllelePresenceGene(String gene) {
    return ALLELE_PRESENCE_GENES.contains(gene);
  }

  /**
   * Returns true if the gene uses activity score to assign phenotype and match to recommendations.
   *
   * @param gene the gene symbol
   * @return true if this gene is an activity score gene
   */
  public static boolean isActivityScoreGene(String gene) {
    return gene != null && ACTIVITY_SCORE_GENES.contains(gene.toUpperCase());
  }

  public static boolean isXChromo(String gene) {
    return StringUtils.isNotBlank(gene) && X_CHROMOSOME_GENES.contains(gene);
  }
}
