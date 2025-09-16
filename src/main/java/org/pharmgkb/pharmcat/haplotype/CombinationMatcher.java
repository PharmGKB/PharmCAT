package org.pharmgkb.pharmcat.haplotype;

import java.util.*;
import com.google.common.base.Splitter;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;


/**
 * This class helps with combination matches.
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
    for (String seq : matchData.getPermutations()) {
      // generate allele map
      SortedMap<Long, String> alleleMap = new TreeMap<>();
      VariantLocus[] refVariants = matchData.getPositions();
      SortedSet<Long> varPositions = new TreeSet<>();
      for (int x = 0; x < refVariants.length; x += 1) {
        long pos = refVariants[x].getPosition();
        String allele = matchData.getAllele(seq, x);
        alleleMap.put(pos, allele);
        if (!refVariants[x].getRef().equals(allele)) {
          varPositions.add(pos);
        }
      }

      // get all possible haplotype matches
      SortedSet<NamedAllele> coveredHaps = new TreeSet<>();
      for (NamedAllele hap : matchData.getHaplotypes()) {
        if (!hap.isReference() && sampleHasNamedAllele(alleleMap, hap)) {
          coveredHaps.add(hap);
        }
      }

      Map<Long, String> partialNames = calculatePartialNames(alleleMap, varPositions, coveredHaps);

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

        if (partialNames.isEmpty()) {
          HaplotypeMatch simpleMatch = new HaplotypeMatch(hap);
          simpleMatch.addSequence(seq);
          matches.add(simpleMatch);
        } else {
          matches.add(new CombinationMatch(refVariants, seq, List.of(hap), partialNames));
        }

      } else {
        List<SortedSet<NamedAllele>> combos = computeViableCombinations(coveredHaps);
        for (SortedSet<NamedAllele> combo : combos) {
          Map<Long, String> partials = partialNames;
          if (combo.size() == coveredHaps.size()) {
            partials = calculatePartialNames(alleleMap, varPositions, combo);
          }
          matches.add(new CombinationMatch(refVariants, seq, combo, partials));
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
  private boolean sampleHasNamedAllele(Map<Long, String> alleleMap, NamedAllele namedAllele) {

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
}
