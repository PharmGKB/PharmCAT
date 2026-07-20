package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.concurrent.TimeUnit;
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
import org.pharmgkb.pharmcat.TestVcfBuilder;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.benchmark.BenchmarkVcfBuilder;
import org.pharmgkb.pharmcat.haplotype.benchmark.PositionsIndex;
import org.pharmgkb.pharmcat.haplotype.model.Result;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * JMH microbenchmark measuring the full per-sample matching path end to end (MatchData construction,
 * marshallHaplotypes, generateSamplePermutations, comparePermutations, and DiplotypeMatcher pairing/scoring), not
 * just index build.
 *
 * <p>Run with: {@code ./gradlew jmh}</p>
 *
 * <p>Each {@link NamedAlleleMatcher} is created once in {@link #setup} and reused across invocations, so after
 * warmup the per-gene {@link HaplotypeCandidateIndex} is already cached and the benchmark isolates the remaining
 * per-sample cost. {@code cyp2c19} exercises the standard reference-defaulted path; {@code dpyd} exercises the
 * lowest-function fallback ladder, including {@link DpydHapB3Matcher}.</p>
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 2)
@Measurement(iterations = 5, time = 2)
@Fork(2)
public class FullCallPathBench {

  private NamedAlleleMatcher m_cyp2c19Matcher;
  private VcfFile m_cyp2c19VcfFile;

  private NamedAlleleMatcher m_dpydMatcher;
  private VcfFile m_dpydVcfFile;

  /** Temporary directory cleaned up after the trial. */
  private Path m_tempDir;


  @Setup(Level.Trial)
  public void setup() throws Exception {
    m_tempDir = Files.createTempDirectory("jmh-fullcall-");
    PositionsIndex positionsIndex = PositionsIndex.getInstance();

    m_cyp2c19VcfFile = buildSampleVcf(positionsIndex, "CYP2C19", "cyp2c19_het.vcf");
    m_cyp2c19Matcher = new NamedAlleleMatcher(TestVcfBuilder.DEFAULT_TEST_ENV,
        new DefinitionReader(List.of(DataManager.getDefinitionFilePath("CYP2C19")), DataManager.DEFAULT_EXEMPTIONS_FILE),
        /*findCombinations=*/false, /*topCandidateOnly=*/true, /*callCyp2d6=*/false);

    m_dpydVcfFile = buildSampleVcf(positionsIndex, "DPYD", "dpyd_het.vcf");
    m_dpydMatcher = new NamedAlleleMatcher(TestVcfBuilder.DEFAULT_TEST_ENV,
        new DefinitionReader(List.of(DataManager.getDefinitionFilePath("DPYD")), DataManager.DEFAULT_EXEMPTIONS_FILE),
        /*findCombinations=*/false, /*topCandidateOnly=*/true, /*callCyp2d6=*/false);
  }


  /**
   * Builds a sample VCF with the first few positions set heterozygous (unphased), forcing a non-trivial diplotype
   * call instead of a trivial all-reference match.
   */
  private VcfFile buildSampleVcf(PositionsIndex positionsIndex, String gene, String fileName) throws Exception {
    DefinitionReader definitionReader = new DefinitionReader(
        List.of(DataManager.getDefinitionFilePath(gene)),
        DataManager.DEFAULT_EXEMPTIONS_FILE);
    BenchmarkVcfBuilder builder = new BenchmarkVcfBuilder(gene, positionsIndex, definitionReader);
    int n = Math.min(6, builder.size());
    for (int i = 0; i < n; i++) {
      builder.set(i, "0/1");
    }
    Path vcfPath = m_tempDir.resolve(fileName);
    builder.write(vcfPath);
    return new VcfFile(vcfPath);
  }


  @TearDown(Level.Trial)
  public void tearDown() throws Exception {
    if (m_tempDir != null) {
      Files.deleteIfExists(m_tempDir.resolve("cyp2c19_het.vcf"));
      Files.deleteIfExists(m_tempDir.resolve("dpyd_het.vcf"));
      Files.deleteIfExists(m_tempDir);
    }
  }


  /**
   * Standard gene, reference-defaulted exact-matching path.
   */
  @Benchmark
  public void cyp2c19(Blackhole bh) throws Exception {
    Result result = m_cyp2c19Matcher.call(m_cyp2c19VcfFile, null, null);
    bh.consume(result);
  }


  /**
   * Lowest-function gene fallback ladder, including the {@link DpydHapB3Matcher} stage.
   */
  @Benchmark
  public void dpyd(Blackhole bh) throws Exception {
    Result result = m_dpydMatcher.call(m_dpydVcfFile, null, null);
    bh.consume(result);
  }
}
