package org.pharmgkb.pharmcat.haplotype;

import java.util.*;
import com.google.common.base.Splitter;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


/**
 * Matches a sample strand to non-overlapping combinations of named alleles and optional partial alleles.
 *
 * <p>This is the research/lowest-function path, not the standard exact matcher. A small core-position index removes
 * impossible named alleles before {@link #sampleHasNamedAllele(Map, NamedAllele)} performs the complete check.</p>
 *
 * @author Mark Woon
 */
public class CombinationMatcher {
  public static final String COMBINATION_JOINER = " + ";
  public static final String COMBINATION_JOINER_REGEX = " \\+ ";
  private static final Splitter COMBINATION_NAME_SPLITTER = Splitter.on(COMBINATION_JOINER).trimResults();
  private final DefinitionFile m_definitionFile;
  private final boolean m_findPartials;


  public static boolean isCombinationName(String name) {
    return name.startsWith("[") && name.contains(COMBINATION_JOINER) && name.endsWith("]");
  }

  public static String extractCombinationName(String name) {
    return name.substring(1, name.length() - 1);
  }

  /**
   * Splits combination name into components.
   * Assumes {@code name} is in combination format (i.e. enclosed within square brackets []).
   */
  public static List<String> splitCombinationName(String name) {
    return COMBINATION_NAME_SPLITTER.splitToList(extractCombinationName(name));
  }



  public CombinationMatcher(DefinitionFile definitionFile, boolean findPartials) {
    m_definitionFile = definitionFile;
    m_findPartials = findPartials;
  }


  /**
   * Compute matches.
   */
  public SortedSet<BaseMatch> compute(MatchData matchData) {

    SortedSet<BaseMatch> matches = new TreeSet<>();
    CandidateIndex candidateIndex = new CandidateIndex(matchData.getHaplotypes());
    for (SamplePermutation permutation : matchData.getPermutations()) {
      // generate allele map
      SortedMap<Long, @Nullable String> alleleMap = new TreeMap<>();
      VariantLocus[] refVariants = matchData.getPositions();
      SortedSet<Long> varPositions = new TreeSet<>();
      for (VariantLocus refVariant : refVariants) {
        long pos = refVariant.getPosition();
        String allele = matchData.getAllele(permutation, pos);
        alleleMap.put(pos, allele);
        if (!refVariant.getRef().equals(allele)) {
          varPositions.add(pos);
        }
      }

      // The index is only a gate. Always run the full core-position check on every candidate it returns.
      SortedSet<NamedAllele> coveredHaps = new TreeSet<>();
      for (NamedAllele hap : candidateIndex.getCandidates(alleleMap)) {
        if (sampleHasNamedAllele(alleleMap, hap)) {
          coveredHaps.add(hap);
        }
      }

      if (coveredHaps.size() <= 1) {
        NamedAllele hap;
        if (coveredHaps.isEmpty()) {
          hap = matchData.getHaplotypes().stream()
              .filter(NamedAllele::isReference)
              .findFirst()
              .orElseThrow();
        } else {
          hap = coveredHaps.first();
        }

        Map<Long, String> partialNames = calculatePartialNames(alleleMap, varPositions, coveredHaps);
        if (partialNames.isEmpty()) {
          HaplotypeMatch simpleMatch = new HaplotypeMatch(hap);
          simpleMatch.addSequence(permutation);
          matches.add(simpleMatch);
        } else {
          matches.add(new CombinationMatch(refVariants, permutation.getSequence(), List.of(hap), partialNames));
        }

      } else {
        List<SortedSet<NamedAllele>> combos = computeViableCombinations(coveredHaps);
        for (SortedSet<NamedAllele> combo : combos) {
          Map<Long, String> partialNames = calculatePartialNames(alleleMap, varPositions, combo);
          matches.add(new CombinationMatch(refVariants, permutation.getSequence(), combo, partialNames));
        }
      }
    }
    return matches;
  }


  private static final Comparator<NamedAllele> sf_numCorePosComparator = Comparator
      .comparingInt((NamedAllele na) -> na.getCorePositions().size())
      .reversed()
      .thenComparing(NamedAllele::getName);

  /**
   * Compute viable combinations of {@link NamedAllele}s.
   * This will take overlapping named alleles and missing positions into account.
   */
  private List<SortedSet<NamedAllele>> computeViableCombinations(Collection<NamedAllele> coveredHaps) {
    SortedSet<NamedAllele> sortedHaps = new TreeSet<>(sf_numCorePosComparator);
    sortedHaps.addAll(coveredHaps);

    List<SortedSet<NamedAllele>> combos = new ArrayList<>();
    for (NamedAllele allele : sortedHaps) {
      boolean added = false;
      int overlapCount = 0;
      for (Set<NamedAllele> combo : combos) {
        boolean overlaps = false;
        for (NamedAllele existingAllele : combo) {
          if (overlaps(existingAllele.getCorePositions(), allele.getCorePositions())) {
            overlaps = true;
            // if both are the same size, it's probably because of a missing position, so don't count it as an overlap
            if (existingAllele.getCorePositions().size() != allele.getCorePositions().size()) {
              overlapCount += 1;
              break;
            }
          }
        }
        if (!overlaps) {
          combo.add(allele);
          added = true;
        }
      }
      if (!added && overlapCount == 0) {
        SortedSet<NamedAllele> newCombo = new TreeSet<>();
        newCombo.add(allele);
        combos.add(newCombo);
      }
    }
    return combos;
  }

  private boolean overlaps(SortedSet<Long> existingPositions, SortedSet<Long> newPositions) {
    for (Long pos : newPositions) {
      if (existingPositions.contains(pos)) {
        return true;
      }
    }
    return false;
  }


  private Map<Long, String> calculatePartialNames(SortedMap<Long, String> alleleMap, SortedSet<Long> varPositions,
      Set<NamedAllele> coveredHaps) {
    if (!m_findPartials) {
      return Collections.emptyMap();
    }
    Map<Long, String> partialNames = new HashMap<>();
    SortedSet<Long> partialPositions = new TreeSet<>(varPositions);
    coveredHaps.forEach(hap -> partialPositions.removeAll(hap.getCorePositions()));
    for (Long pos : partialPositions) {
      partialNames.put(pos, m_definitionFile.getVariantForPosition(pos).getHgvsForVcfAllele(alleleMap.get(pos)));
    }
    return partialNames;
  }


  /**
   * Checks if a sample has all the alleles for the specified {@code namedAllele}.
   */
  private boolean sampleHasNamedAllele(Map<Long, @Nullable String> alleleMap, NamedAllele namedAllele) {

    for (long pos : namedAllele.getCorePositions()) {
      String sampleAllele = alleleMap.get(pos);
      if (sampleAllele == null) {
        return false;
      }
      String expectedAllele = Objects.requireNonNull(namedAllele.getAllele(pos));
      if (namedAllele.isWobble(pos)) {
        if (!Iupac.lookup(expectedAllele).getBases().contains(sampleAllele)) {
          return false;
        }
      } else if (!expectedAllele.equals(sampleAllele)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Gates each non-reference haplotype on its first core position.
   *
   * <p>The gate may return false positives because it checks only one position; that is intentional and is corrected
   * by {@link CombinationMatcher#sampleHasNamedAllele(Map, NamedAllele)}. It must not return false negatives: exact
   * and wobble gates therefore use the same observed position and allele semantics as the full check.</p>
   */
  private static class CandidateIndex {
    private final Map<Long, Map<String, SortedSet<NamedAllele>>> m_exactCandidates = new HashMap<>();
    private final Map<Long, List<WobbleCandidate>> m_wobbleCandidates = new HashMap<>();
    private final SortedSet<NamedAllele> m_alwaysCheckCandidates = new TreeSet<>();

    CandidateIndex(SortedSet<NamedAllele> haplotypes) {
      for (NamedAllele hap : haplotypes) {
        if (hap.isReference()) {
          continue;
        }
        if (hap.getCorePositions().isEmpty()) {
          m_alwaysCheckCandidates.add(hap);
          continue;
        }
        long position = hap.getCorePositions().first();
        String allele = Objects.requireNonNull(hap.getAllele(position));
        if (hap.isWobble(position)) {
          m_wobbleCandidates.computeIfAbsent(position, p -> new ArrayList<>())
              .add(new WobbleCandidate(hap, allele));
        } else {
          m_exactCandidates.computeIfAbsent(position, p -> new HashMap<>())
              .computeIfAbsent(allele, a -> new TreeSet<>())
              .add(hap);
        }
      }
    }

    SortedSet<NamedAllele> getCandidates(Map<Long, @Nullable String> alleleMap) {
      SortedSet<NamedAllele> candidates = new TreeSet<>(m_alwaysCheckCandidates);
      for (Long position : alleleMap.keySet()) {
        String allele = alleleMap.get(position);
        Map<String, SortedSet<NamedAllele>> exactCandidates = m_exactCandidates.get(position);
        if (exactCandidates != null) {
          SortedSet<NamedAllele> exactMatches = exactCandidates.get(allele);
          if (exactMatches != null) {
            candidates.addAll(exactMatches);
          }
        }
        List<WobbleCandidate> wobbleCandidates = m_wobbleCandidates.get(position);
        if (wobbleCandidates != null && allele != null) {
          for (WobbleCandidate candidate : wobbleCandidates) {
            if (candidate.matches(allele)) {
              candidates.add(candidate.haplotype());
            }
          }
        }
      }
      return candidates;
    }
  }

  private record WobbleCandidate(NamedAllele haplotype, String allele) {

    boolean matches(String sampleAllele) {
      return Iupac.lookup(allele).getBases().contains(sampleAllele);
    }
  }
}
