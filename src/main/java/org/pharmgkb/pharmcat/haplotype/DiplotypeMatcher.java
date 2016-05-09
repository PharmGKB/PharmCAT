package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
import javax.annotation.Nonnull;
import com.google.common.collect.Lists;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


/**
 * This is the main class responsible for calling diplotypes.
 *
 * @author Mark Woon
 */
public class DiplotypeMatcher {
  private MatchData m_dataset;


  public DiplotypeMatcher(@Nonnull MatchData dataset) {

    m_dataset = dataset;
  }


  public List<DiplotypeMatch> compute() {

    // compare sample permutations to haplotypes
    SortedSet<HaplotypeMatch> matches = comparePermutations();

    if (m_dataset.permutations.size() == 1 && matches.size() == 1) {
      // sample is homozygous for all positions and it matches a single allele,
      // so we need to return that as a diplotype
      HaplotypeMatch hm = matches.first();
      DiplotypeMatch dm = new DiplotypeMatch(hm, hm, m_dataset);
      String seq = m_dataset.permutations.iterator().next();
      dm.addSequencePair(new String[] { seq, seq });
      return Lists.newArrayList(dm);
    }

    // find matched pairs
    return determinePairs(matches);
  }


  /**
   * Compares a sample's allele permutations to haplotype definitions and return matches.
   */
  protected @Nonnull SortedSet<HaplotypeMatch> comparePermutations() {

    Set<HaplotypeMatch> haplotypeMatches = m_dataset.haplotypes.stream()
        .map(HaplotypeMatch::new)
        .collect(Collectors.toSet());

    for (String p : m_dataset.permutations) {
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
   * @param haplotypeMatches the matches that were found via {@link #comparePermutations()}
   */
  protected @Nonnull List<DiplotypeMatch> determinePairs(@Nonnull SortedSet<HaplotypeMatch> haplotypeMatches) {

    SortedMap<NamedAllele, HaplotypeMatch> hapMap = new TreeMap<>();
    for (HaplotypeMatch hm : haplotypeMatches) {
      hapMap.put(hm.getHaplotype(), hm);
    }

    // possible pairs from what got matched
    List<List<NamedAllele>> pairs = CombinationUtil.generatePerfectPairs(hapMap.keySet());

    List<DiplotypeMatch> matches = new ArrayList<>();
    for (List<NamedAllele> pair : pairs) {
      NamedAllele hap1 = pair.get(0);
      HaplotypeMatch hm1 = hapMap.get(hap1);
      if (hm1 == null) {
        continue;
      }
      NamedAllele hap2 = pair.get(1);
      HaplotypeMatch hm2 = hapMap.get(hap2);
      if (hm2 == null) {
        continue;
      }

      if (hap1 == hap2 && hm1.getSequences().size() == 1) {
        // cannot call homozygous unless more than one sequence matches
        continue;
      }

      Set<String[]> sequencePairs = findSequencePairs(hm1, hm2);
      if (!sequencePairs.isEmpty()) {
        DiplotypeMatch dm = new DiplotypeMatch(hm1, hm2, m_dataset);
        sequencePairs.stream().forEach(dm::addSequencePair);
        matches.add(dm);
      }
    }
    Collections.sort(matches);
    return matches;
  }


  /**
   * Finds valid complementary pairs of sample's alleles for possible diplotype match.
   */
  private Set<String[]> findSequencePairs(@Nonnull HaplotypeMatch hm1, @Nonnull HaplotypeMatch hm2) {

    Set<String[]> sequencePairs = new HashSet<>();
    for (String seq1 : hm1.getSequences()) {
      for (String seq2 : hm2.getSequences()) {
        if (isViableComplement(seq1, seq2)) {
          sequencePairs.add(new String[] { seq1, seq2 });
        }
      }
    }
    return sequencePairs;
  }


  /**
   * Checks whether the two sequences is complementary based on sample alleles.
   */
  private boolean isViableComplement(@Nonnull String sequence1, @Nonnull String sequence2) {

    String[] seq1 = sequence1.split(";");
    String[] seq2 = sequence2.split(";");

    for (int x = 0; x < seq1.length; x += 1) {
      String[] s1 = seq1[x].split(":");
      String[] s2 = seq2[x].split(":");
      SampleAllele sampleAllele = m_dataset.geneSampleMap.get(Integer.valueOf(s1[0]));
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
}
