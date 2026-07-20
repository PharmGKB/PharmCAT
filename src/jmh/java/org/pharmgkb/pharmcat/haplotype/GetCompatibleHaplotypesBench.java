package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.concurrent.TimeUnit;
import org.jspecify.annotations.Nullable;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Fork;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.TearDown;
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.benchmark.BenchmarkVcfBuilder;
import org.pharmgkb.pharmcat.haplotype.benchmark.PositionsIndex;
import org.pharmgkb.pharmcat.reporter.TextConstants;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * JMH microbenchmark comparing the current per-position full-scan {@code getCompatibleHaplotypes}
 * build strategy against a proposed sparse-index build, on a live all-reference RYR1 {@link MatchData}.
 *
 * <p>Run with: {@code ./gradlew jmh}</p>
 *
 * <p>Both benchmarks build the complete per-position BitSet index from scratch each invocation
 * (no cross-invocation state reuse), so the comparison reflects real cold-build cost.</p>
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 2)
@Measurement(iterations = 5, time = 2)
@Fork(2)
public class GetCompatibleHaplotypesBench {

  /** Number of permutation positions (= number of positions available for the sample). */
  private int m_numPositions;

  /** Number of callable haplotypes. */
  private int m_numHaplotypes;

  /**
   * Expected alleles in [haplotypeIndex][positionIndex] order, post-defaulting.
   * Mirrors {@code m_haplotypeAlleles.get(h)[p]} after {@code initializeCandidateIndex} runs.
   */
  private @Nullable String[][] m_haplotypeAlleles;

  /**
   * Raw (pre-defaulting) alleles in [haplotypeIndex][positionIndex] order.
   * Used by sparseImpl to expose null (wildcard) slots.
   */
  private @Nullable String[][] m_rawHaplotypeAlleles;

  /** The callable haplotype list in stable index order. */
  private List<NamedAllele> m_haplotypeIndex;

  /** Observed alleles from the first SamplePermutation, in permutation-position order. */
  private @Nullable String[] m_observedAlleles;

  /** Temporary directory cleaned up after the trial. */
  private Path m_tempDir;


  @Setup(Level.Trial)
  public void setup() throws Exception {
    m_tempDir = Files.createTempDirectory("jmh-gch-");

    // Build an all-reference RYR1 VCF (0/0 everywhere — default, no edits needed).
    String gene = "RYR1";
    DefinitionReader definitionReader = new DefinitionReader(
        List.of(DataManager.getDefinitionFilePath(gene)),
        DataManager.DEFAULT_EXEMPTIONS_FILE);

    PositionsIndex positionsIndex = PositionsIndex.getInstance();
    BenchmarkVcfBuilder builder = new BenchmarkVcfBuilder(gene, positionsIndex, definitionReader);
    // No set() calls — every position defaults to 0/0 (all-reference).
    Path vcfPath = m_tempDir.resolve("ryr1_allref.vcf");
    builder.write(vcfPath);

    // Read the VCF exactly as NamedAlleleMatcher.call() does.
    VcfFile vcfFile = new VcfFile(vcfPath);
    // assumeReference=true, findCombinations=false (the callAssumingReference path).
    VcfReader vcfReader = vcfFile.getReader(definitionReader, null, /*findCombinations=*/false);
    SortedMap<String, SampleAllele> alleleMap = vcfReader.getAlleleMap();
    String sampleId = vcfReader.getSampleId();

    // Build MatchData via the same sequence as NamedAlleleMatcher.prepareMatchData(assumeReference=true, findCombinations=false).
    VariantLocus[] allPositions = definitionReader.getPositions(gene);
    MatchData matchData = new MatchData(sampleId, gene, alleleMap, allPositions,
        /*extraPositions=*/null, /*exemption=*/null);
    matchData.marshallHaplotypes(gene, definitionReader.getHaplotypes(gene), /*findCombinations=*/false,
        /*assumeReference=*/true);
    matchData.generateSamplePermutations();

    // Extract the permutation positions (numeric/sorted order — same as m_permutationPositions).
    VariantLocus[] permutationPositions = Arrays.stream(matchData.getPositions())
        .sorted()
        .toArray(VariantLocus[]::new);
    m_numPositions = permutationPositions.length;

    // Extract first permutation's observed alleles (in permutation-position order).
    SamplePermutation firstPerm = matchData.getPermutations().iterator().next();
    m_observedAlleles = firstPerm.getAllelesForMatching();

    // Build the stable haplotype index (mirrors initializeCandidateIndex's m_haplotypeIndex).
    m_haplotypeIndex = new ArrayList<>(matchData.getHaplotypes());
    m_numHaplotypes = m_haplotypeIndex.size();

    // Build reference alleles array for defaulting (mirrors getAllelesForMatching logic).
    NamedAllele referenceHaplotype = matchData.getHaplotypes().stream()
        .filter(NamedAllele::isReference)
        .findAny()
        .orElseThrow(() -> new IllegalStateException("RYR1 has no reference haplotype"));
    @Nullable String[] referenceAlleles = referenceHaplotype.getAlleles(permutationPositions);

    // Precompute haplotypeAlleles[h][p] — post-defaulting, mirrors m_haplotypeAlleles.
    m_haplotypeAlleles = new String[m_numHaplotypes][];
    m_rawHaplotypeAlleles = new String[m_numHaplotypes][];
    for (int h = 0; h < m_numHaplotypes; h++) {
      NamedAllele hap = m_haplotypeIndex.get(h);
      @Nullable String[] raw = hap.getAlleles(permutationPositions);
      m_rawHaplotypeAlleles[h] = raw;
      if (hap.isReference()) {
        // Reference haplotype: no defaulting applied (getAllelesForMatching returns raw as-is).
        m_haplotypeAlleles[h] = raw;
      } else {
        // Non-reference: null slots defaulted to reference allele (mirrors applyReferenceDefaults).
        @Nullable String[] defaulted = new String[m_numPositions];
        for (int p = 0; p < m_numPositions; p++) {
          if (raw[p] == null) {
            String refAllele = referenceAlleles[p];
            if (Iupac.isWobble(refAllele)) {
              defaulted[p] = permutationPositions[p].getRef();
            } else {
              defaulted[p] = refAllele;
            }
          } else {
            defaulted[p] = raw[p];
          }
        }
        m_haplotypeAlleles[h] = defaulted;
      }
    }

    // Correctness check: verify that both impls produce identical BitSets for every position.
    verifyEquivalence();
  }


  @TearDown(Level.Trial)
  public void tearDown() throws Exception {
    if (m_tempDir != null) {
      Path vcf = m_tempDir.resolve("ryr1_allref.vcf");
      Files.deleteIfExists(vcf);
      Files.deleteIfExists(m_tempDir);
    }
  }


  /**
   * Assert that currentImpl and sparseImpl produce identical BitSets at every position.
   * Runs once per Trial so we catch any wobble-expansion bug before benchmarking.
   */
  private void verifyEquivalence() {
    BitSet[] currentResults = buildCurrentIndex();
    BitSet[] sparseResults = buildSparseIndex();
    for (int p = 0; p < m_numPositions; p++) {
      if (!currentResults[p].equals(sparseResults[p])) {
        throw new IllegalStateException(
            "current[" + p + "] != sparse[" + p + "] for observed=" + m_observedAlleles[p]
                + ": current=" + currentResults[p] + " sparse=" + sparseResults[p]);
      }
    }
  }


  /**
   * Current implementation: for each position, walk all haplotypes and call matchesAllele.
   * Mirrors what {@code MatchData.getCompatibleHaplotypes} does on a cold cache (first permutation pass).
   */
  @Benchmark
  public void currentImpl(Blackhole bh) {
    BitSet[] results = buildCurrentIndex();
    for (BitSet bs : results) {
      bh.consume(bs.cardinality());
    }
  }


  /**
   * Proposed sparse implementation: build per-position wildcard + specific-allele BitSets,
   * then union them for the observed allele at each position.
   */
  @Benchmark
  public void sparseImpl(Blackhole bh) {
    BitSet[] results = buildSparseIndex();
    for (BitSet bs : results) {
      bh.consume(bs.cardinality());
    }
  }


  /**
   * Builds one result BitSet per position using the current full-scan strategy.
   * Each result bit p is the set of haplotypes compatible with observedAlleles[p].
   */
  private BitSet[] buildCurrentIndex() {
    BitSet[] results = new BitSet[m_numPositions];
    for (int p = 0; p < m_numPositions; p++) {
      BitSet compatible = new BitSet(m_numHaplotypes);
      String observed = m_observedAlleles[p];
      for (int h = 0; h < m_numHaplotypes; h++) {
        if (m_haplotypeIndex.get(h).matchesAllele(m_haplotypeAlleles[h][p], observed)) {
          compatible.set(h);
        }
      }
      results[p] = compatible;
    }
    return results;
  }


  /**
   * Builds one result BitSet per position using the proposed sparse-index strategy.
   *
   * <p>Build phase: iterate haplotypes × defining positions only (raw != null).
   * For each (h, p) where raw is non-null, expand wobble semantics exactly as
   * {@code NamedAllele.matchesAllele} does, and add h to S(p, each expanded base).
   * Wildcards (raw == null) go into W(p).</p>
   *
   * <p>Query phase: for each position p, result = clone(W(p)).or(S(p, observed[p])).</p>
   */
  @SuppressWarnings("unchecked")
  private BitSet[] buildSparseIndex() {
    // W(p): haplotypes with null raw allele at position p (match anything).
    BitSet[] wildcards = new BitSet[m_numPositions];
    // S(p, allele): haplotypes with a specific raw allele (after wobble expansion) at position p.
    Map<String, BitSet>[] specifics = new HashMap[m_numPositions];
    for (int p = 0; p < m_numPositions; p++) {
      wildcards[p] = new BitSet(m_numHaplotypes);
      specifics[p] = new HashMap<>();
    }

    // Build phase: walk haplotypes × defining positions.
    for (int h = 0; h < m_numHaplotypes; h++) {
      @Nullable String[] rawAlleles = m_rawHaplotypeAlleles[h];
      for (int p = 0; p < m_numPositions; p++) {
        String raw = rawAlleles[p];
        if (raw == null) {
          wildcards[p].set(h);
        } else {
          // Expand wobble semantics to match NamedAllele.matchesAllele exactly.
          for (String base : expandAllele(raw)) {
            specifics[p].computeIfAbsent(base, k -> new BitSet(m_numHaplotypes)).set(h);
          }
        }
      }
    }

    // Query phase: union wildcard set with the specific set for the observed allele.
    BitSet[] results = new BitSet[m_numPositions];
    for (int p = 0; p < m_numPositions; p++) {
      BitSet result = (BitSet) wildcards[p].clone();
      BitSet specific = specifics[p].get(m_observedAlleles[p]);
      if (specific != null) {
        result.or(specific);
      }
      results[p] = result;
    }
    return results;
  }


  /**
   * Expands a raw haplotype allele into the set of observed-allele strings it matches,
   * mirroring {@code NamedAllele.matchesAllele} semantics exactly.
   *
   * <ul>
   *   <li>REPEAT_WOBBLE_DELIMITER present → split tokens (each token is one matchable string).</li>
   *   <li>length==1 and IUPAC ambiguity → {@code Iupac.getBases()} (e.g. "R" → ["A","G"]).</li>
   *   <li>length==1 and not ambiguity → {@code Iupac.getRegex()} (e.g. "A" → "A", "-" → "del").</li>
   *   <li>Otherwise (multi-char, no wobble delimiter) → raw itself.</li>
   * </ul>
   */
  private static List<String> expandAllele(String raw) {
    if (raw.contains(TextConstants.REPEAT_WOBBLE_DELIMITER)) {
      return Arrays.asList(raw.split(TextConstants.REPEAT_WOBBLE_DELIMITER));
    }
    if (raw.length() == 1) {
      Iupac iupac = Iupac.lookup(raw);
      if (iupac.isAmbiguity()) {
        return iupac.getBases();
      }
      return List.of(iupac.getRegex());
    }
    return List.of(raw);
  }
}
