package org.cpic.haplotype;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;


/**
 * @author Mark Woon
 */
public class ComparisonUtil {


  public static SortedSet<HaplotypeMatch> comparePermutations(Set<String> permutations, List<Haplotype> haplotypes) {

    Set<HaplotypeMatch> haplotypeMatches = haplotypes.stream()
        .map(HaplotypeMatch::new)
        .collect(Collectors.toSet());

    for (String p : permutations) {
      for (HaplotypeMatch hm : haplotypeMatches) {
        hm.match(p);
      }
    }

    return haplotypeMatches.stream()
        .filter(h -> !h.getSequences().isEmpty())
        .collect(Collectors.toCollection(TreeSet::new));
  }



  public static void determinePairs(Set<HaplotypeMatch> haplotypeMatches, List<List<Haplotype>> pairs) {

    Map<Haplotype, HaplotypeMatch> hapMap = new HashMap<>();
    for (HaplotypeMatch hm : haplotypeMatches) {
      hapMap.put(hm.getHaplotype(), hm);
    }

    for (List<Haplotype> pair : pairs) {
      Haplotype hap1 = pair.get(0);
      HaplotypeMatch hm1 = hapMap.get(hap1);
      if (hm1 == null) {
        continue;
      }
      Haplotype hap2 = pair.get(1);
      HaplotypeMatch hm2 = hapMap.get(hap2);
      if (hm2 == null) {
        continue;
      }

      if (hap1 == hap2 && hm1.getSequences().size() == 1) {
        // cannot call homozygous unless more than one sequence matches
        continue;
      }

      System.out.println(pair);
      System.out.println(hap1);
      System.out.println(hm1.getSequences());
      if (hap1 != hap2) {
        System.out.println(hap2);
        System.out.println(hm2.getSequences());
      }
      System.out.println();
      System.out.println("------------");
      System.out.println();
    }
  }
}
