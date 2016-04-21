package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


/**
 * @author Mark Woon
 */
public class ComparisonUtil {


  /**
   * Compares a sample's allele permutations to haplotype definitions and return matches.
   */
  public static @Nonnull SortedSet<HaplotypeMatch> comparePermutations(@Nonnull Set<String> permutations,
      @Nonnull Collection<Haplotype> haplotypes) {

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
   *
   * @param haplotypeMatches the matches that were found via {@link #comparePermutations(Set, Collection)}
   */
  public static List<DiplotypeMatch> determinePairs(@Nonnull SortedMap<Integer, SampleAllele> haplotypeSampleMap,
      @Nonnull SortedSet<HaplotypeMatch> haplotypeMatches) {

    SortedMap<Haplotype, HaplotypeMatch> hapMap = new TreeMap<>();
    for (HaplotypeMatch hm : haplotypeMatches) {
      hapMap.put(hm.getHaplotype(), hm);
    }

    // possible pairs from what got matched
    List<List<Haplotype>> pairs = CombinationUtil.generatePerfectPairs(hapMap.keySet());

    List<DiplotypeMatch> matches = new ArrayList<>();
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

      Set<String[]> sequencePairs = findSequencePairs(haplotypeSampleMap, hm1, hm2);
      if (!sequencePairs.isEmpty()) {
        DiplotypeMatch dm = new DiplotypeMatch(hm1, hm2);
        sequencePairs.stream().forEach(dm::addSequencePair);
        matches.add(dm);
      }
    }
    return matches;
  }


  private static Set<String[]> findSequencePairs(@Nonnull SortedMap<Integer, SampleAllele> haplotypeSampleMap,
      @Nonnull HaplotypeMatch hm1, @Nonnull HaplotypeMatch hm2) {

    Set<String[]> sequencePairs = new HashSet<>();
    for (String seq1 : hm1.getSequences()) {
      for (String seq2 : hm2.getSequences()) {
        if (isViableComplement(haplotypeSampleMap, seq1, seq2)) {
          sequencePairs.add(new String[] { seq1, seq2 });
        }
      }
    }
    return sequencePairs;
  }


  /**
   * Checks whether the two sequences is complementary based on sample alleles.
   */
  private static boolean isViableComplement(@Nonnull SortedMap<Integer, SampleAllele> haplotypeSampleMap,
      @Nonnull String sequence1, @Nonnull String sequence2) {

    String[] seq1 = sequence1.split(";");
    String[] seq2 = sequence2.split(";");

    for (int x = 0; x < seq1.length; x += 1) {
      String[] s1 = seq1[x].split(":");
      String[] s2 = seq2[x].split(":");
      SampleAllele sampleAllele = haplotypeSampleMap.get(Integer.valueOf(s1[0]));
      if (sampleAllele == null) {
        throw new IllegalStateException("Missing sample for haplotype position " + s1[0]);
      }
      if (sampleAllele.getAllele1().equals(sampleAllele.getAllele2())) {
        // expecting homozygous
        if (!s1[1].equals(s2[1])) {
          return false;
        }
      } else {
        // expecting heterozygous
        if (s1[1].equals(s2[1])) {
          return false;
        }
      }
    }

    return true;
  }


  public static void printMatchPairs(List<DiplotypeMatch> matches) {

    for (DiplotypeMatch pair : matches) {
      HaplotypeMatch hm1 = pair.getHaplotype1();
      HaplotypeMatch hm2 = pair.getHaplotype2();

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
