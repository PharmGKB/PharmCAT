package org.cpic.haplotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import javax.annotation.Nonnull;
import com.google.common.collect.Sets;


/**
 * @author Mark Woon
 */
public class CombinationUtil {

  /**
   * Builds permutations for given alleles based on phasing.
   */
  public static Set<String> generatePermutations(@Nonnull List<SampleAllele> alleles) {
    return generatePermutations(alleles, 0, true, "");
  }


  /**
   * Builds permutations for given variants based on phasing.
   */
  private static Set<String> generatePermutations(@Nonnull List<SampleAllele> sampleAlleles, int position, boolean isFirst,
      @Nonnull String alleleSoFar) {

    if (position >= sampleAlleles.size()) {
      return Sets.newHashSet(alleleSoFar);
    }
    SampleAllele allele = sampleAlleles.get(position);

    Set<String> alleles = new HashSet<>();
    if (allele.isPhased()) {
      generatePermutations(sampleAlleles, position + 1, isFirst, appendAllele(alleleSoFar, allele, isFirst));
    } else {
      alleles.addAll(generatePermutations(sampleAlleles, position + 1, isFirst, appendAllele(alleleSoFar, allele, true)));
      alleles.addAll(generatePermutations(sampleAlleles, position + 1, isFirst, appendAllele(alleleSoFar, allele, false)));
    }
    return alleles;
  }

  private static String appendAllele(String alleleSoFar, SampleAllele allele, boolean isFirst) {
    StringBuilder sb = new StringBuilder()
        .append(alleleSoFar)
        .append(allele.getPosition())
        .append(":");
    if (isFirst) {
      sb.append(allele.getAllele1());
    } else {
      sb.append(allele.getAllele2());
    }
    sb.append(";");
    return sb.toString();
  }


  public static <T> List<List<T>> generatePerfectPairs(Collection<T> data) {

    List<T> list;
    if (data instanceof List) {
      //noinspection unchecked
      list = (List)data;
    } else {
      list = new ArrayList<>(data);
    }
    List<List<T>> rez = new ArrayList<>();
    for (int x = 0; x < list.size(); x++) {
      for (int y = x; y < list.size(); y++) {
        rez.add(Arrays.asList(list.get(x), list.get(y)));
      }
    }
    return rez;
  }
}
