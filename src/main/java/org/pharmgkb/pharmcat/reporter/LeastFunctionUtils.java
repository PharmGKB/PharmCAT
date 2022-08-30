package org.pharmgkb.pharmcat.reporter;

import java.util.Set;
import com.google.common.collect.ImmutableSet;
import org.pharmgkb.pharmcat.haplotype.model.GeneCall;


/**
 * This class centralizes least-function support functions.
 *
 * @author Mark Woon
 */
public class LeastFunctionUtils {
  private static final Set<String> LEAST_FUNCTION = ImmutableSet.of("DPYD");

  /**
   * Gets whether this gene should use the "least-function" algorithm for determining alleles.
   */
  public static boolean useLeastFunction(String gene) {
    return LEAST_FUNCTION.contains(gene);
  }


  /**
   * Checks if a least function gene has a {@link GeneCall} with a true diplotype.
   */
  public static boolean hasTrueDiplotype(GeneCall geneCall) {
    if (!geneCall.isEffectivelyPhased()) {
      return false;
    }
    if (geneCall.getDiplotypes().size() > 1) {
      throw new IllegalStateException("Least function gene cannot have more than 1 diplotype");
    }
    return geneCall.getDiplotypes().size() == 1;
  }
}
