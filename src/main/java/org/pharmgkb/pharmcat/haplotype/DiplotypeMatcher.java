package org.pharmgkb.pharmcat.haplotype;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;
import java.util.stream.Collectors;
import com.google.common.collect.Sets;
import com.google.common.collect.SortedSetMultimap;
import com.google.common.collect.TreeMultimap;
import org.jspecify.annotations.Nullable;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.model.BaseMatch;
import org.pharmgkb.pharmcat.haplotype.model.CombinationMatch;
import org.pharmgkb.pharmcat.haplotype.model.DiplotypeMatch;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.util.HaplotypeNameComparator;


/**
 * Pairs single-strand matches into diplotypes that can explain the sample genotype.
 *
 * <p>Single-strand matching is delegated to {@link MatchData} or {@link CombinationMatcher}. This class validates
 * complementary strand pairs, scores exact named-allele diplotypes, and applies top-candidate filtering. Structured
 * {@link SamplePermutation} metadata is preferred, with encoded sequence parsing retained for legacy and synthetic
 * combination matches.</p>
 *
 * @author Mark Woon
 */
public class DiplotypeMatcher {
  private final MatchData m_dataset;
  private final DefinitionFile m_definitionFile;
  private final boolean m_unphasedPriorityMode;
  private final VariantLocus[] m_positions;
  /** VCF positions aligned with m_positions for structured complement checks. */
  private final long[] m_vcfPositions;
  /** Sample zygosity aligned with m_positions and m_vcfPositions. */
  private final boolean[] m_isHomozygous;
  private final boolean m_timing;


  public DiplotypeMatcher(Env env, MatchData dataset) {
    this(env, dataset, false);
  }

  public DiplotypeMatcher(Env env, MatchData dataset, boolean timing) {
    m_dataset = dataset;
    m_timing = timing;
    m_positions = dataset.getPositions();
    m_vcfPositions = new long[m_positions.length];
    m_isHomozygous = new boolean[m_positions.length];
    for (int x = 0; x < m_positions.length; x += 1) {
      long position = m_positions[x].getPosition();
      m_vcfPositions[x] = position;
      m_isHomozygous[x] = dataset.getSampleAllele(position).isHomozygous();
    }
    m_definitionFile = env.getDefinitionReader().getDefinitionFile(dataset.getGene());
    DefinitionExemption exemption = env.getDefinitionReader().getExemption(dataset.getGene());
    m_unphasedPriorityMode = !dataset.isEffectivelyPhased() &&
        exemption != null && exemption.hasUnphasedDiplotypePriorities();
  }


  /**
   * Computes diplotypes, enabling partial-allele matching whenever combination matching is enabled.
   */
  public SortedSet<DiplotypeMatch> compute(boolean findCombinations, boolean topCandidateOnly) {
    return compute(findCombinations, findCombinations, topCandidateOnly);
  }

  public SortedSet<DiplotypeMatch> compute(boolean findCombinations, boolean findPartials, boolean topCandidateOnly) {
    if (findCombinations && topCandidateOnly) {
      throw new IllegalStateException("Cannot get top candidate only when using combinations!");
    }
    String context = m_dataset.getGene() + " DiplotypeMatcher.compute(findCombinations=" + findCombinations + ")";
    long totalStart = MatcherTimings.start(m_timing);

    SortedSet<BaseMatch> matches;
    if (findCombinations) {
      long stageStart = MatcherTimings.start(m_timing);
      matches = new CombinationMatcher(m_definitionFile, findPartials)
          .compute(m_dataset);
      MatcherTimings.print(m_timing, context, "CombinationMatcher.compute", stageStart);

    } else {
      // compare sample permutations to haplotypes
      long stageStart = MatcherTimings.start(m_timing);
      List<HaplotypeMatch> haplotypeMatches = new ArrayList<>(m_dataset.comparePermutations());
      MatcherTimings.print(m_timing, context, "MatchData.comparePermutations", stageStart);
      if (haplotypeMatches.isEmpty()) {
        MatcherTimings.print(m_timing, context, "total", totalStart);
        return Collections.emptySortedSet();
      }
      matches = new TreeSet<>(haplotypeMatches);
    }

    List<DiplotypeMatch> pairs;
    long stageStart = MatcherTimings.start(m_timing);
    if (m_dataset.getPermutationCount() == 1) {
      pairs = determineHomozygousPairs(matches);
    } else {
      // find matched pairs
      pairs = determineHeterozygousPairs(matches, findCombinations);
    }
    MatcherTimings.print(m_timing, context, "determine diplotype pairs", stageStart);

    if (!findCombinations) {
      stageStart = MatcherTimings.start(m_timing);
      for (DiplotypeMatch dm : pairs) {
        // score is based on the best scoring pair of sequences for this diplotype
        int highestScore = 0;
        for (String[] seqPair : dm.getSequences()) {
          int score = scoreForSequencePair(dm, seqPair);
          if (score > highestScore) {
            highestScore = score;
          }
        }
        dm.setScore(highestScore);
      }
      MatcherTimings.print(m_timing, context, "score diplotype pairs", stageStart);
    }

    // TODO(markwoon): if combinations, and phased, and we have more than one match, it's probably because of wobbles
    // TODO(markwoon): if there are wobbles, should we use the top candidate only?

    // using Collections.sort() throws an exception, so use SortedSet instead
    SortedSet<DiplotypeMatch> sortedPairs = new TreeSet<>(pairs);
    if (topCandidateOnly && !m_unphasedPriorityMode && sortedPairs.size() > 1) {
      int topScore = sortedPairs.first().getScore();
      SortedSet<DiplotypeMatch> topMatches = sortedPairs.stream()
          .filter(dm -> dm.getScore() == topScore)
          .collect(Collectors.toCollection(TreeSet::new));
      MatcherTimings.print(m_timing, context, "total", totalStart);
      return topMatches;
    }
    MatcherTimings.print(m_timing, context, "total", totalStart);
    return sortedPairs;
  }

  private int scoreForSequencePair(DiplotypeMatch dm, String[] seqPair) {
    BaseMatch m1 = dm.getHaplotype1();
    BaseMatch m2 = dm.getHaplotype2();
    boolean isHomozygous = false;
    int m2Score = 0;
    if (m2 != null) {
      if (m1.getName().equals(m2.getName())) {
        isHomozygous = true;
        m2Score = scoreForBaseMatch(m2, new String[] { seqPair[1] });
      } else {
        m2Score = scoreForBaseMatch(m2, seqPair);
      }
    }
    if (isHomozygous) {
      return scoreForBaseMatch(m1, new String[] {seqPair[0] }) + m2Score;
    } else {
      return scoreForBaseMatch(m1, seqPair) + m2Score;
    }
  }

  private int scoreForBaseMatch(BaseMatch hapMatch, String[] seqPair) {
    Set<String> sequences = sequencesForBaseMatch(hapMatch, seqPair);
    List<SamplePermutation> permutations = new ArrayList<>(sequences.size());
    for (String sequence : sequences) {
      SamplePermutation permutation = hapMatch.getSequencePermutation(sequence);
      if (permutation == null) {
        return hapMatch.getHaplotype().scoreForSample(m_dataset, sequences);
      }
      permutations.add(permutation);
    }
    return scoreForSample(hapMatch.getHaplotype(), permutations);
  }

  private int scoreForSample(NamedAllele haplotype, List<SamplePermutation> permutations) {
    if (haplotype.getWobblePositions().isEmpty()) {
      return haplotype.getScore();
    }
    int score = haplotype.getScore();
    for (Long position : haplotype.getWobblePositions()) {
      VariantLocus vl = getPosition(position);
      int numRefs = 0;
      for (SamplePermutation permutation : permutations) {
        String allele = m_dataset.getAllele(permutation, position);
        if (Objects.equals(allele, vl.getRef())) {
          numRefs += 1;
        }
      }
      // if all alleles at position are ref, don't score this wobble
      if (numRefs == permutations.size()) {
        score -= 1;
      }
    }
    return score;
  }

  private VariantLocus getPosition(long position) {
    for (VariantLocus vl : m_positions) {
      if (vl.getPosition() == position) {
        return vl;
      }
    }
    throw new IllegalArgumentException("Unknown position: " + position);
  }

  private Set<String> sequencesForBaseMatch(BaseMatch hapMatch, String[] seqPair) {
    Set<String> seqs = hapMatch.getSequences();
    if (seqs.size() == 1) {
      return seqs;
    }
    seqs = Sets.intersection(seqs, new HashSet<>(Arrays.asList(seqPair)));
    if (seqs.isEmpty()) {
      throw new IllegalStateException("NamedAllele sequences do not match Diplotype sequence!");
    }
    return seqs;
  }


  /**
   * Determine possible diplotypes given a set of {@link HaplotypeMatch}'s when sample is homozygous at all positions.
   *
   * @param haplotypeMatches the matches that were found via {@link MatchData#comparePermutations()}
   */
  private List<DiplotypeMatch> determineHomozygousPairs(SortedSet<BaseMatch> haplotypeMatches) {

    String seq = m_dataset.getPermutations().iterator().next().getSequence();
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

    // map haplotype name to HaplotypeMatches (i.e., sequences)
    SortedSetMultimap<String, BaseMatch> hapMap = TreeMultimap.create();
    for (BaseMatch hm : haplotypeMatches) {
      hapMap.put(hm.getName(), hm);
    }

    // possible pairs from what got matched
    List<String> sortedNames = new ArrayList<>(hapMap.keySet());
    sortedNames.sort(HaplotypeNameComparator.getComparator());
    List<List<String>> pairs = CombinationUtil.generatePerfectPairs(sortedNames);

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
        SamplePermutation permutation1 = hm1.getSequencePermutation(seq1);
        SamplePermutation permutation2 = hm2.getSequencePermutation(seq2);
        boolean viable;
        // Combination/partial and reconstructed DPYD matches may not retain structured permutation metadata.
        if (permutation1 != null && permutation2 != null) {
          viable = isViableComplement(permutation1, permutation2);
        } else {
          viable = isViableComplement(seq1, seq2);
        }
        if (viable) {
          sequencePairs.add(new String[] { seq1, seq2 });
        }
      }
    }
    return sequencePairs;
  }


  /**
   * Checks whether two encoded sequences are complementary based on sample alleles.
   * This is the compatibility fallback for matches without structured permutation metadata.
   */
  private boolean isViableComplement(String sequence1, String sequence2) {

    for (int x = 0; x < m_positions.length; x += 1) {
      String a1 = m_dataset.getAllele(sequence1, x);
      String a2 = m_dataset.getAllele(sequence2, x);
      if (!matchesSampleZygosity(a1, a2, m_isHomozygous[x])) {
        return false;
      }
    }

    return true;
  }

  /**
   * Checks whether two structured permutations are complementary based on sample alleles.
   */
  private boolean isViableComplement(SamplePermutation permutation1, SamplePermutation permutation2) {

    for (int x = 0; x < m_positions.length; x += 1) {
      long position = m_vcfPositions[x];
      String a1 = m_dataset.getAllele(permutation1, position);
      String a2 = m_dataset.getAllele(permutation2, position);
      if (!matchesSampleZygosity(a1, a2, m_isHomozygous[x])) {
        return false;
      }
    }

    return true;
  }


  /**
   * Applies the complement rule at one position. Homozygous samples require equal strand alleles; heterozygous
   * samples require different strand alleles. {@link Objects#equals(Object, Object)} keeps both paths null-safe.
   */
  private static boolean matchesSampleZygosity(@Nullable String allele1, @Nullable String allele2,
      boolean sampleIsHomozygous) {
    return sampleIsHomozygous == Objects.equals(allele1, allele2);
  }
}
