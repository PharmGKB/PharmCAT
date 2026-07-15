package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import com.google.common.base.Preconditions;
import org.jspecify.annotations.Nullable;


/**
 * @author Mark Woon
 */
public class CombinationUtil {

  /**
   * Builds permutations for given alleles based on phasing.
   */
  public static Set<String> generatePermutations(List<SampleAllele> alleles) {
    Set<String> rez = new HashSet<>();
    for (SamplePermutation permutation : generatePermutationData(alleles)) {
      rez.add(permutation.getSequence());
    }
    return rez;
  }


  /**
   * Builds internal allele-array permutations for given alleles based on phasing.
   */
  static Set<SamplePermutation> generatePermutationData(List<SampleAllele> alleles) {
    Preconditions.checkNotNull(alleles);
    Preconditions.checkArgument(!alleles.isEmpty(), "No alleles to generate permutations for");

    boolean isS1Blank = true;
    boolean isS2Blank = true;
    boolean hasPhaseSets = false;
    for (SampleAllele sa : alleles) {
      if (sa.getPhaseSet() != null) {
        hasPhaseSets = true;
      }
      if (sa.getAllele1() != null) {
        isS1Blank = false;
      }
      if (sa.getAllele2() != null) {
        isS2Blank = false;
      }
      if (!isS1Blank && !isS2Blank && hasPhaseSets) {
        break;
      }
    }
    long[] positions = alleles.stream()
        .mapToLong(SampleAllele::getPosition)
        .toArray();
    Set<SamplePermutation> rez = new HashSet<>();
    String[] alleleSoFar = new String[alleles.size()];
    if (!isS1Blank) {
      Map<Integer, Boolean> phaseSets = hasPhaseSets ? new HashMap<>() : null;
      generatePermutations(alleles, 0, isS2Blank, true, alleleSoFar, positions, phaseSets, rez);
    }
    if (!isS2Blank) {
      Map<Integer, Boolean> phaseSets = hasPhaseSets ? new HashMap<>() : null;
      generatePermutations(alleles, 0, isS1Blank, false, alleleSoFar, positions, phaseSets, rez);
    }
    if (rez.isEmpty()) {
      throw new IllegalStateException("No permutations generated from " + alleles.size() + " alleles");
    }
    return rez;
  }


  /**
   * Builds permutations for given variants based on phasing.
   */
  private static void generatePermutations(List<SampleAllele> sampleAlleles, int position, boolean isHaploid,
      boolean firstAllele, @Nullable String[] alleleSoFar, long[] positions,
      @Nullable Map<Integer, Boolean> phaseSets, Set<SamplePermutation> permutations) {

    if (position >= sampleAlleles.size()) {
      permutations.add(new SamplePermutation(alleleSoFar, positions));
      return;
    }
    SampleAllele allele = sampleAlleles.get(position);

    if (allele.isEffectivelyPhased() || isHaploid) {
      alleleSoFar[position] = firstAllele ? allele.getComputedAllele1() : allele.getComputedAllele2();
      generatePermutations(sampleAlleles, position + 1, isHaploid, firstAllele, alleleSoFar, positions,
          phaseSets, permutations);
    } else if (allele.getPhaseSet() != null) {
      //noinspection DataFlowIssue
      if (phaseSets.containsKey(allele.getPhaseSet())) {
        if (phaseSets.get(allele.getPhaseSet())) {
          alleleSoFar[position] = allele.getComputedAllele1();
          generatePermutations(sampleAlleles, position + 1, false, firstAllele, alleleSoFar, positions,
              phaseSets, permutations);
        } else {
          alleleSoFar[position] = allele.getComputedAllele2();
          generatePermutations(sampleAlleles, position + 1, false, firstAllele, alleleSoFar, positions,
              phaseSets, permutations);
        }
      } else {
        // initial PS key
        // in phase set
        Map<Integer, Boolean> ps1 = new HashMap<>(phaseSets);
        ps1.put(allele.getPhaseSet(), true);
        alleleSoFar[position] = allele.getComputedAllele1();
        generatePermutations(sampleAlleles, position + 1, false, firstAllele, alleleSoFar, positions,
            ps1, permutations);
        // out of phase set
        Map<Integer, Boolean> ps2 = new HashMap<>(phaseSets);
        ps2.put(allele.getPhaseSet(), false);
        alleleSoFar[position] = allele.getComputedAllele2();
        generatePermutations(sampleAlleles, position + 1, false, firstAllele, alleleSoFar, positions,
            ps2, permutations);
      }
    } else {
      alleleSoFar[position] = allele.getComputedAllele1();
      generatePermutations(sampleAlleles, position + 1, false, firstAllele, alleleSoFar, positions,
          phaseSets, permutations);
      alleleSoFar[position] = allele.getComputedAllele2();
      generatePermutations(sampleAlleles, position + 1, false, firstAllele, alleleSoFar, positions,
          phaseSets, permutations);
    }
  }


  public static <T> List<List<T>> generatePerfectPairs(Collection<T> data) {

    List<T> list;
    if (data instanceof List) {
      list = (List<T>)data;
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
