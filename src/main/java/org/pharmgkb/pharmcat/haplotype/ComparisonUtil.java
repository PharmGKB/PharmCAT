package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.collect.Lists;


/**
 * @author Mark Woon
 */
public class ComparisonUtil {


  /**
   * Compares a sample's allele permutations to haplotype definitions and return matches.
   */
  public static SortedSet<HaplotypeMatch> comparePermutations(Set<String> permutations,
      Collection<Haplotype> haplotypes) {

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


  /**
   * Determine possible diplotypes given a set of {@link HaplotypeMatch}'s.
   */
  public static List<List<HaplotypeMatch>> determinePairs(Set<HaplotypeMatch> haplotypeMatches,
      List<List<Haplotype>> pairs) {

    Map<Haplotype, HaplotypeMatch> hapMap = new HashMap<>();
    for (HaplotypeMatch hm : haplotypeMatches) {
      hapMap.put(hm.getHaplotype(), hm);
    }

    List<List<HaplotypeMatch>> matches = new ArrayList<>();
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
      matches.add(Lists.newArrayList(hm1, hm2));
    }
    return matches;
  }


  public static void printMatchPairs(List<List<HaplotypeMatch>> matches) {

    for (List<HaplotypeMatch> pair : matches) {
      HaplotypeMatch hm1 = pair.get(0);
      HaplotypeMatch hm2 = pair.get(1);

      System.out.println(pair);
      System.out.println(hm1.getHaplotype());
      System.out.println(hm1.getSequences());
      if (hm1.getHaplotype() != hm2.getHaplotype()) {
        System.out.println(hm2.getHaplotype());
        System.out.println(hm2.getSequences());
      }
      System.out.println();
      System.out.println("------------");
      System.out.println();
    }
  }
}
