package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
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
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.benchmark.BenchmarkVcfBuilder;
import org.pharmgkb.pharmcat.haplotype.benchmark.PositionsIndex;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * JMH microbenchmark comparing the current array-based {@code isViableComplement} implementation
 * against the old HashMap-based implementation, isolated from end-to-end matcher noise.
 *
 * <p>Run with: {@code ./gradlew jmh}</p>
 *
 * <p>Both benchmarks set up the same {@link MatchData} and {@link SamplePermutation} pair that
 * the {@code combinations/CYP2D6/het=6} scenario produces, then call the inner loop repeatedly.</p>
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.NANOSECONDS)
@Warmup(iterations = 5, time = 2)
@Measurement(iterations = 5, time = 2)
@Fork(2)
public class IsViableComplementBench {

  /** Live MatchData built from the combinations/CYP2D6/het=6 scenario. */
  private MatchData m_matchData;

  /** First permutation for the complement check. */
  private SamplePermutation m_perm1;

  /** Second permutation (different from m_perm1). */
  private SamplePermutation m_perm2;

  /**
   * Precomputed permutation-position index array, mirroring
   * {@code DiplotypeMatcher.m_permIndexByPositionOrder}.
   */
  private int[] m_permIndexByPositionOrder;

  /**
   * Sample zygosity per position in m_positions order, mirroring
   * {@code DiplotypeMatcher.m_isHomozygous}.
   */
  private boolean[] m_isHomozygous;

  /** Positions in definition/output order, mirroring {@code DiplotypeMatcher.m_positions}. */
  private VariantLocus[] m_positions;

  /** Temporary directory cleaned up after the trial. */
  private Path m_tempDir;


  @Setup(Level.Trial)
  public void setup() throws Exception {
    m_tempDir = Files.createTempDirectory("jmh-isviable-");

    // Build the combinations/CYP2D6/het=6 VCF using the same builder the Tier-1 harness uses.
    String gene = "CYP2D6";
    DefinitionReader definitionReader = new DefinitionReader(
        List.of(DataManager.getDefinitionFilePath(gene)),
        DataManager.DEFAULT_EXEMPTIONS_FILE);

    PositionsIndex positionsIndex = PositionsIndex.getInstance();
    BenchmarkVcfBuilder builder = new BenchmarkVcfBuilder(gene, positionsIndex, definitionReader);
    // set first 6 positions heterozygous (unphased) — same as Scenarios.combinations("CYP2D6", 6)
    int n = Math.min(6, builder.size());
    for (int i = 0; i < n; i++) {
      builder.set(i, "0/1");
    }
    Path vcfPath = m_tempDir.resolve("cyp2d6_het6.vcf");
    builder.write(vcfPath);

    // Read the VCF and build an alleleMap, exactly as NamedAlleleMatcher.call() does.
    VcfFile vcfFile = new VcfFile(vcfPath);
    // findCombinations=true to match the CombinationMatcher code path
    VcfReader vcfReader = vcfFile.getReader(definitionReader, null, /*findCombinations=*/true);
    SortedMap<String, SampleAllele> alleleMap = vcfReader.getAlleleMap();
    String sampleId = vcfReader.getSampleId();

    // Build MatchData via the same path as NamedAlleleMatcher.initializeCallData(…, false, true).
    VariantLocus[] allPositions = definitionReader.getPositions(gene);
    m_matchData = new MatchData(sampleId, gene, alleleMap, allPositions,
        /*extraPositions=*/null, /*exemption=*/null);
    m_matchData.marshallHaplotypes(gene, definitionReader.getHaplotypes(gene), /*findCombinations=*/true);
    // No defaultMissingAllelesToReference for combinations mode (assumeReference=false).
    m_matchData.generateSamplePermutations();

    // Extract two DIFFERENT permutations.
    List<SamplePermutation> perms = new ArrayList<>(m_matchData.getPermutations());
    if (perms.size() < 2) {
      throw new IllegalStateException(
          "Expected >= 2 permutations for het scenario, got " + perms.size());
    }
    m_perm1 = perms.get(0);
    m_perm2 = perms.get(1);

    // Mirror DiplotypeMatcher constructor: precompute index array and zygosity array.
    m_positions = m_matchData.getPositions();
    m_permIndexByPositionOrder = new int[m_positions.length];
    m_isHomozygous = new boolean[m_positions.length];
    for (int x = 0; x < m_positions.length; x++) {
      long position = m_positions[x].getPosition();
      m_permIndexByPositionOrder[x] = m_matchData.getPermutationIndex(position);
      m_isHomozygous[x] = m_matchData.getSampleAllele(position).isHomozygous();
    }
  }


  @TearDown(Level.Trial)
  public void tearDown() throws Exception {
    // Delete the temp VCF file we created.
    if (m_tempDir != null) {
      Path vcf = m_tempDir.resolve("cyp2d6_het6.vcf");
      Files.deleteIfExists(vcf);
      Files.deleteIfExists(m_tempDir);
    }
  }


  /**
   * Current (array-based) implementation: precomputed {@code int[]} index lookup.
   * Mirrors {@code DiplotypeMatcher.isViableComplement(SamplePermutation, SamplePermutation)}.
   */
  @Benchmark
  public boolean arrayBased() {
    @Nullable String[] alleles1 = m_perm1.getAllelesForMatching();
    @Nullable String[] alleles2 = m_perm2.getAllelesForMatching();

    for (int x = 0; x < m_positions.length; x++) {
      int permIdx = m_permIndexByPositionOrder[x];
      String a1 = alleles1[permIdx];
      String a2 = alleles2[permIdx];
      if (!matchesSampleZygosity(a1, a2, m_isHomozygous[x])) {
        return false;
      }
    }
    return true;
  }


  /**
   * Old (HashMap-based) implementation: {@code m_permutationPositionIndex.get(pos)} per position.
   * Mirrors the pre-refactor {@code DiplotypeMatcher.isViableComplement(SamplePermutation, SamplePermutation)}.
   */
  @Benchmark
  public boolean hashMapBased() {
    for (int x = 0; x < m_positions.length; x++) {
      long position = m_positions[x].getPosition();
      // MatchData.getAllele(SamplePermutation, long) does m_permutationPositionIndex.get(position)
      String a1 = m_matchData.getAllele(m_perm1, position);
      String a2 = m_matchData.getAllele(m_perm2, position);
      if (!matchesSampleZygosity(a1, a2, m_isHomozygous[x])) {
        return false;
      }
    }
    return true;
  }


  private static boolean matchesSampleZygosity(@Nullable String allele1, @Nullable String allele2,
      boolean sampleIsHomozygous) {
    return sampleIsHomozygous == Objects.equals(allele1, allele2);
  }
}
