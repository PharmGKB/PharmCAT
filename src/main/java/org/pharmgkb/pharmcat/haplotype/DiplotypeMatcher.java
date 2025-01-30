package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.checkerframework.checker.nullness.qual.Nullable;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


/**
 * This is the main class responsible for calling diplotypes.
 *
 * @author Mark Woon
 */
public class DiplotypeMatcher {
  private final MatchData m_dataset;


  public DiplotypeMatcher(MatchData dataset) {
    m_dataset = dataset;
  }


  public SortedSet<DiplotypeMatch> compute(boolean findCombinations) {
    return compute(findCombinations, false);
  }

  public SortedSet<DiplotypeMatch> compute(boolean findCombinations, boolean topCandidateOnly) {
    return compute(findCombinations, findCombinations, topCandidateOnly, false);
  }

  /**
   *
   * @param boostComboScores false by default, which means scoring will prefer fewer combos; if true, scoring will
   * prefer more combos
   */
  public SortedSet<DiplotypeMatch> compute(boolean findCombinations, boolean findPartials, boolean topCandidateOnly,
      boolean boostComboScores) {

    // compare sample permutations to haplotypes
    List<HaplotypeMatch> haplotypeMatches = new ArrayList<>(m_dataset.comparePermutations());
    if (haplotypeMatches.isEmpty() && !findCombinations) {
      return Collections.emptySortedSet();
    }

    SortedSet<BaseMatch> matches = new TreeSet<>();
    if (findCombinations) {
      List<CombinationMatch> combinationMatches = new ArrayList<>();
      for (int x = 0; x < haplotypeMatches.size(); x += 1) {
        HaplotypeMatch hm = haplotypeMatches.get(x);
        for (String seq : hm.getSequences()) {
          // to support partial matching, each HaplotypeMatch can only have 1 sequence
          HaplotypeMatch individualHapMatch = new HaplotypeMatch(hm.getHaplotype());
          individualHapMatch.addSequence(seq);
          individualHapMatch.finalizeCombinationHaplotype(m_dataset, findPartials);
          matches.add(individualHapMatch);

          CombinationMatch combinationMatch = new CombinationMatch(m_dataset.getPositions(), hm.getHaplotype(), seq);
          combinationMatches.addAll(calculateCombinations(haplotypeMatches, x + 1, combinationMatch));
        }
      }

      // finalize the combination match and check for partials (if findPartials = true)
      for (CombinationMatch combinationMatch : combinationMatches) {
        combinationMatch.finalizeCombinationHaplotype(m_dataset, findPartials);
        matches.add(combinationMatch);
      }
      if (findPartials) {
        // add remaining permutations as off-reference partial match
        NamedAllele referenceHaplotype = m_dataset.getHaplotypes().stream()
            .filter(NamedAllele::isReference)
            .findFirst()
            .orElseThrow(() -> new IllegalStateException("No reference haplotype!"));
        Set<String> matchedSequences = matches.stream()
            .flatMap(bm -> bm.getSequences().stream())
            .collect(Collectors.toSet());
        for (String seq : m_dataset.getPermutations()) {
          if (matchedSequences.contains(seq)) {
            continue;
          }
          matches.add(new CombinationMatch(m_dataset, referenceHaplotype, seq));
        }
      }

    } else {
      matches.addAll(haplotypeMatches);
    }

    List<DiplotypeMatch> pairs;
    if (m_dataset.getPermutations().size() == 1) {
      pairs = determineHomozygousPairs(matches);
    } else {
      // find matched pairs
      pairs = determineHeterozygousPairs(matches, findCombinations);
    }

    // final scoring
    if (findCombinations) {
      // weed out shorter pairs
      int maxCombos = 0;
      for (DiplotypeMatch dm : pairs) {
        maxCombos = maxComboBonus(dm.getHaplotype1(), maxCombos);
        maxCombos = maxComboBonus(dm.getHaplotype2(), maxCombos);
      }
      for (DiplotypeMatch dm : pairs) {
        dm.setScore(newComboScore(dm.getHaplotype1(), maxCombos, boostComboScores) +
            newComboScore(dm.getHaplotype2(), maxCombos, boostComboScores));
      }
    } else {
      for (DiplotypeMatch dm : pairs) {
        BaseMatch m1 = dm.getHaplotype1();
        BaseMatch m2 = dm.getHaplotype2();
        int m2Score = 0;
        if (m2 != null) {
          m2Score = m2.getHaplotype().scoreForSample(m_dataset, m2.getSequences());
        }
        dm.setScore(m1.getHaplotype().scoreForSample(m_dataset, m1.getSequences()) + m2Score);
      }
    }

    // using Collections.sort() throws an exception, so use SortedSet instead
    SortedSet<DiplotypeMatch> sortedPairs = new TreeSet<>(pairs);
    if (topCandidateOnly && sortedPairs.size() > 1) {
      int topScore = sortedPairs.first().getScore();
      return sortedPairs.stream()
          .filter(dm -> dm.getScore() == topScore)
          .collect(Collectors.toCollection(TreeSet::new));
    }
    return sortedPairs;
  }

  private int maxComboBonus(@Nullable BaseMatch match, int curMax) {
    if (match == null) {
      return curMax;
    }
    if (match instanceof CombinationMatch cm) {
      if (cm.getNumCombinations() > curMax) {
        return cm.getNumCombinations();
      }
    }
    return curMax;
  }

  private int newComboScore(@Nullable BaseMatch match, int maxBonus, boolean boostCombos) {
    if (match == null) {
      return 0;
    }
    NamedAllele na = match.getHaplotype();
    int score = na.scoreForSample(m_dataset, match.getSequences()) - na.getNumPartials();
    int bonus;
    if (boostCombos) {
      // prefer more combos
      bonus = match.getHaplotype().getNumCombinations();
    } else {
      // default behavior - prefer fewer combos
      bonus = maxBonus - match.getHaplotype().getNumCombinations();
    }
    return score + bonus;
  }


  private List<CombinationMatch> calculateCombinations(List<HaplotypeMatch> matches, int position,
      CombinationMatch combinationMatch) {

    List<CombinationMatch> combinationMatches = new ArrayList<>();
    for (int x = position; x < matches.size(); x += 1) {
      HaplotypeMatch haplotypeMatch = matches.get(x);
      for (String seq : haplotypeMatch.getSequences()) {
        if (combinationMatch.canMerge(haplotypeMatch.getHaplotype(), seq)) {
          CombinationMatch newCombo = new CombinationMatch(combinationMatch);
          newCombo.merge(haplotypeMatch.getHaplotype());
          combinationMatches.add(newCombo);
          combinationMatches.addAll(calculateCombinations(matches, x + 1, newCombo));
        }
      }
    }
    return combinationMatches;
  }


  /**
   * Determine possible diplotypes given a set of {@link HaplotypeMatch}'s when sample is homozygous at all positions.
   *
   * @param haplotypeMatches the matches that were found via {@link MatchData#comparePermutations()}
   */
  private List<DiplotypeMatch> determineHomozygousPairs(SortedSet<BaseMatch> haplotypeMatches) {

    String seq = m_dataset.getPermutations().iterator().next();
    List<DiplotypeMatch> matches = new ArrayList<>();
    if (haplotypeMatches.size() == 1) {
      // matched a single haplotype: need to return that as either homozygous diplotype or haploid
      BaseMatch hm1 = haplotypeMatches.first();
      BaseMatch hm2 = null;
      String[] sequencePair;
      if (m_dataset.isHaploid()) {
        sequencePair = new String[] {seq};
      } else {
        hm2 = hm1;
        sequencePair = new String[] {seq, seq};
      }
      DiplotypeMatch dm = new DiplotypeMatch(hm1, hm2, m_dataset);
      dm.addSequencePair(sequencePair);
      matches.add(dm);
    } else {
      // return all possible pairings of matched haplotypes
      List<List<BaseMatch>> pairs = CombinationUtil.generatePerfectPairs(haplotypeMatches);
      for (List<BaseMatch> pair : pairs) {
        DiplotypeMatch dm = new DiplotypeMatch(pair.get(0), pair.get(1), m_dataset);
        dm.addSequencePair(new String[]{ seq, seq });
        matches.add(dm);
      }
    }
    return matches;
  }



  /**
   * Determine possible diplotypes given a set of {@link HaplotypeMatch}'s when sample is heterozygous at (at least) one
   * position.
   *
   * @param haplotypeMatches the matches that were found via {@link MatchData#comparePermutations()}
   */
  private List<DiplotypeMatch> determineHeterozygousPairs(SortedSet<BaseMatch> haplotypeMatches, boolean findCombinations) {

    SortedSetMultimap<String, BaseMatch> hapMap = TreeMultimap.create();
    for (BaseMatch hm : haplotypeMatches) {
      hapMap.put(hm.getName(), hm);
    }

    // possible pairs from what got matched
    List<List<String>> pairs = CombinationUtil.generatePerfectPairs(new TreeSet<>(hapMap.keySet()));

    List<DiplotypeMatch> matches = new ArrayList<>();
    for (List<String> pair : pairs) {
      String name1 = pair.get(0);
      SortedSet<BaseMatch> hm1s = hapMap.get(name1);
      String name2 = pair.get(1);
      SortedSet<BaseMatch> hm2s = hapMap.get(name2);

      if (name1.equals(name2)) {
        // hm1s and hm2s collections are the same
        if (hm1s.size() == 1) {
          // if HaplotypeMatch and only has one sequence, cannot be homozygous
          // if CombinationMatch, it will only have 1 sequence, so cannot be homozygous
          if (hm1s.first().getSequences().size() == 1) {
            continue;
          }
        } else {
          if (findCombinations) {
            if (hm1s.first() instanceof CombinationMatch && hm1s.size() == 2) {
              hm2s = new TreeSet<>();
              hm2s.add(hm1s.first());
              // must create new hm1s to avoid messing with backing hapMap
              SortedSet<BaseMatch> new1 = new TreeSet<>();
              new1.add(hm1s.last());
              hm1s = new1;
            }
          }
        }
      }

      for (BaseMatch m1 : hm1s) {
        for (BaseMatch m2 : hm2s) {
          Set<String[]> sequencePairs = findSequencePairs(m1, m2);
          if (!sequencePairs.isEmpty()) {
            DiplotypeMatch dm = new DiplotypeMatch(m1, m2, m_dataset);
            sequencePairs.forEach(dm::addSequencePair);
            matches.add(dm);
          }
        }
      }
    }
    return matches;
  }


  /**
   * Finds valid complementary pairs of sample's alleles for possible diplotype match.
   */
  private Set<String[]> findSequencePairs(BaseMatch hm1, BaseMatch hm2) {

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
  private boolean isViableComplement(String sequence1, String sequence2) {

    for (int x = 0; x < m_dataset.getPositions().length; x += 1) {
      String a1 = m_dataset.getAllele(sequence1, x);
      String a2 = m_dataset.getAllele(sequence2, x);
      SampleAllele sampleAllele = m_dataset.getSampleAllele(m_dataset.getPositions()[x].getPosition());
      if (sampleAllele.isHomozygous()) {
        // expecting homozygous
        if (!a1.equals(a2)) {
          return false;
        }
      } else {
        // expecting heterozygous
        if (a1.equals(a2)) {
          return false;
        }
      }
    }

    return true;
  }
}
