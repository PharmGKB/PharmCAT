package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;
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


  public List<DiplotypeMatch> compute(boolean findCombinations) {
    return compute(findCombinations, false);
  }

  public List<DiplotypeMatch> compute(boolean findCombinations, boolean topCandidateOnly) {

    // compare sample permutations to haplotypes
    List<HaplotypeMatch> haplotypeMatches = new ArrayList<>(m_dataset.comparePermutations());
    if (haplotypeMatches.size() == 0 && !findCombinations) {
      return Collections.emptyList();
    }

    SortedSet<BaseMatch> matches = new TreeSet<>();
    if (findCombinations) {
      List<CombinationMatch> combinationMatches = new ArrayList<>();
      for (int x = 0; x < haplotypeMatches.size(); x += 1) {
        HaplotypeMatch hm = haplotypeMatches.get(x);
        for (String seq : hm.getSequences()) {
          // to support partial matching each HaplotypeMatch can only have 1 sequence
          HaplotypeMatch individualHapMatch = new HaplotypeMatch(hm.getHaplotype());
          individualHapMatch.addSequence(seq);
          individualHapMatch.finalizeCombinationHaplotype(m_dataset);
          matches.add(individualHapMatch);

          CombinationMatch combinationMatch = new CombinationMatch(m_dataset.getPositions(), hm.getHaplotype(), seq);
          combinationMatches.addAll(calculateCombinations(haplotypeMatches, x + 1, combinationMatch));
        }
      }

      // check for partials
      for (CombinationMatch combinationMatch : combinationMatches) {
        combinationMatch.finalizeCombinationHaplotype(m_dataset);
        matches.add(combinationMatch);
      }

      // add off-reference partial match
      // this needs to come after check for partials and added matches do not need to have finalizeHaplotype called
      if (m_dataset.getPermutations().size() <= 2) {
        Optional<NamedAllele> opt = m_dataset.getHaplotypes().stream()
            .filter(NamedAllele::isReference)
            .findFirst();
        if (opt.isPresent()) {
          NamedAllele reference = opt.get();
          for (String seq : m_dataset.getPermutations()) {
            Optional<BaseMatch> match = matches.stream()
                .filter(m -> m.getSequences().contains(seq))
                .findAny();
            if (match.isEmpty()) {
              // only generate if sequence has no other match
              matches.add(new CombinationMatch(m_dataset, reference, seq));
            }
          }
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
      pairs = determineHeterozygousPairs(matches);
    }
    if (findCombinations) {
      // weed out shorter pairs
      int maxCombos = 0;
      for (DiplotypeMatch dm : pairs) {
        maxCombos = maxComboBonus(dm.getHaplotype1(), maxCombos);
        maxCombos = maxComboBonus(dm.getHaplotype2(), maxCombos);
      }
      for (DiplotypeMatch dm : pairs) {
        dm.setScore(newComboScore(dm.getHaplotype1(), maxCombos) +
            newComboScore(dm.getHaplotype2(), maxCombos));
      }
    }
    Collections.sort(pairs);

    if (topCandidateOnly && pairs.size() > 1) {
      int topScore = pairs.get(0).getScore();
      return pairs.stream()
          .filter(dm -> dm.getScore() == topScore)
          .collect(Collectors.toList());
    }
    return pairs;
  }

  private int maxComboBonus(BaseMatch match, int curMax) {
    if (match instanceof CombinationMatch cm) {
      if (cm.getNumCombinations() > curMax) {
        return cm.getNumCombinations();
      }
    }
    return curMax;
  }

  private int newComboScore(BaseMatch match, int maxBonus) {
    NamedAllele na = match.getHaplotype();
    int score = na.getScore() - na.getNumPartials();
    int bonus = maxBonus - match.getHaplotype().getNumCombinations();
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
      // matched a single haplotype: need to return that as a diplotype
      BaseMatch hm = haplotypeMatches.first();
      DiplotypeMatch dm = new DiplotypeMatch(hm, hm, m_dataset);
      dm.addSequencePair(new String[]{ seq, seq });
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
  private List<DiplotypeMatch> determineHeterozygousPairs(SortedSet<BaseMatch> haplotypeMatches) {

    SortedMap<String, BaseMatch> hapMap = new TreeMap<>();
    for (BaseMatch hm : haplotypeMatches) {
      hapMap.put(hm.getName(), hm);
    }

    // possible pairs from what got matched
    List<List<String>> pairs = CombinationUtil.generatePerfectPairs(hapMap.keySet());

    List<DiplotypeMatch> matches = new ArrayList<>();
    for (List<String> pair : pairs) {
      String name1 = pair.get(0);
      BaseMatch hm1 = hapMap.get(name1);

      String name2 = pair.get(1);
      BaseMatch hm2 = hapMap.get(name2);

      if (name1.equals(name2) && hm1.getSequences().size() == 1) {
        // cannot call homozygous unless more than one sequence matches
        continue;
      }

      Set<String[]> sequencePairs = findSequencePairs(hm1, hm2);
      if (!sequencePairs.isEmpty()) {
        DiplotypeMatch dm = new DiplotypeMatch(hm1, hm2, m_dataset);
        sequencePairs.forEach(dm::addSequencePair);
        matches.add(dm);
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
      if (sampleAllele.getAllele1().equals(sampleAllele.getAllele2())) {
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
