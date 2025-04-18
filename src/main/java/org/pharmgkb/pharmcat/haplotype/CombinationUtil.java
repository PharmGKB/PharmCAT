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
import com.google.common.collect.Sets;
import org.checkerframework.checker.nullness.qual.Nullable;


/**
 * @author Mark Woon
 */
public class CombinationUtil {

  /**
   * Builds permutations for given alleles based on phasing.
   */
  public static Set<String> generatePermutations(List<SampleAllele> alleles) {
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
    Set<String> rez = new HashSet<>();
    if (!isS1Blank) {
      Map<Integer, Boolean> phaseSets = hasPhaseSets ? new HashMap<>() : null;
      rez.addAll(generatePermutations(alleles, 0, isS2Blank, true, "", phaseSets));
    }
    if (!isS2Blank) {
      Map<Integer, Boolean> phaseSets = hasPhaseSets ? new HashMap<>() : null;
      rez.addAll(generatePermutations(alleles, 0, isS1Blank, false, "", phaseSets));
    }
    if (rez.isEmpty()) {
      throw new IllegalStateException("No permutations generated from " + alleles.size() + " alleles");
    }
    return rez;
  }


  /**
   * Builds permutations for given variants based on phasing.
   */
  private static Set<String> generatePermutations(List<SampleAllele> sampleAlleles, int position, boolean isHaploid,
      boolean firstAllele, String alleleSoFar, @Nullable Map<Integer, Boolean> phaseSets) {

    if (position >= sampleAlleles.size()) {
      return Sets.newHashSet(alleleSoFar);
    }
    SampleAllele allele = sampleAlleles.get(position);

    Set<String> alleles = new HashSet<>();
    if (allele.isEffectivelyPhased() || isHaploid) {
      alleles.addAll(generatePermutations(sampleAlleles, position + 1, isHaploid, firstAllele,
          appendAllele(alleleSoFar, allele, firstAllele), phaseSets));
    } else if (allele.getPhaseSet() != null) {
      //noinspection DataFlowIssue
      if (phaseSets.containsKey(allele.getPhaseSet())) {
        if (phaseSets.get(allele.getPhaseSet())) {
          alleles.addAll(generatePermutations(sampleAlleles, position + 1, false, firstAllele,
              appendAllele(alleleSoFar, allele, true), phaseSets));
        } else {
          alleles.addAll(generatePermutations(sampleAlleles, position + 1, false, firstAllele,
              appendAllele(alleleSoFar, allele, false), phaseSets));
        }
      } else {
        // initial PS key
        // in phase set
        Map<Integer, Boolean> ps1 = new HashMap<>(phaseSets);
        ps1.put(allele.getPhaseSet(), true);
        alleles.addAll(generatePermutations(sampleAlleles, position + 1, false, firstAllele,
            appendAllele(alleleSoFar, allele, true), ps1));
        // out of phase set
        Map<Integer, Boolean> ps2 = new HashMap<>(phaseSets);
        ps2.put(allele.getPhaseSet(), false);
        alleles.addAll(generatePermutations(sampleAlleles, position + 1, false, firstAllele,
            appendAllele(alleleSoFar, allele, false), ps2));
      }
    } else {
      alleles.addAll(generatePermutations(sampleAlleles, position + 1, false, firstAllele,
          appendAllele(alleleSoFar, allele, true), phaseSets));
      alleles.addAll(generatePermutations(sampleAlleles, position + 1, false, firstAllele,
          appendAllele(alleleSoFar, allele, false), phaseSets));
    }
    return alleles;
  }

  private static String appendAllele(String alleleSoFar, SampleAllele allele, boolean firstAllele) {
    StringBuilder sb = new StringBuilder()
        .append(alleleSoFar)
        .append(allele.getPosition())
        .append(":");
    if (firstAllele) {
      sb.append(allele.getComputedAllele1());
    } else {
      sb.append(allele.getComputedAllele2());
    }
    sb.append(";");
    return sb.toString();
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
