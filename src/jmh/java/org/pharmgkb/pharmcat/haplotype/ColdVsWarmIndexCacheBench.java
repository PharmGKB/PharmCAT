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
 * JMH microbenchmark measuring the payoff of caching the per-gene {@link HaplotypeCandidateIndex} across samples.
 *
 * <p>Run with: {@code ./gradlew jmh}</p>
 *
 * <p>{@code coldCache} constructs a new {@link NamedAlleleMatcher} on every invocation, mirroring today's batch
 * behavior where {@code Pipeline} creates a new matcher per sample and discards its index cache. {@code warmCache}
 * reuses a single matcher instance created once in {@link #setup}, mirroring a hoisted cross-sample cache.</p>
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 2)
@Measurement(iterations = 5, time = 2)
@Fork(2)
public class ColdVsWarmIndexCacheBench {

  /** Gene definitions restricted to a single gene, matching how the matcher is scoped in production. */
  private DefinitionReader m_definitionReader;

  /** All-reference RYR1 VCF, reused for both benchmarks. */
  private VcfFile m_vcfFile;

  /** Matcher created once and reused across invocations (the "warm cache" case). */
  private NamedAlleleMatcher m_warmMatcher;

  /** Temporary directory cleaned up after the trial. */
  private Path m_tempDir;


  @Setup(Level.Trial)
  public void setup() throws Exception {
    m_tempDir = Files.createTempDirectory("jmh-cache-");

    String gene = "RYR1";
    m_definitionReader = new DefinitionReader(
        List.of(DataManager.getDefinitionFilePath(gene)),
        DataManager.DEFAULT_EXEMPTIONS_FILE);

    PositionsIndex positionsIndex = PositionsIndex.getInstance();
    BenchmarkVcfBuilder builder = new BenchmarkVcfBuilder(gene, positionsIndex, m_definitionReader);
    // No set() calls — every position defaults to 0/0 (all-reference).
    Path vcfPath = m_tempDir.resolve("ryr1_allref.vcf");
    builder.write(vcfPath);
    m_vcfFile = new VcfFile(vcfPath);

    m_warmMatcher = new NamedAlleleMatcher(TestVcfBuilder.DEFAULT_TEST_ENV, m_definitionReader,
        /*findCombinations=*/false, /*topCandidateOnly=*/true, /*callCyp2d6=*/false);
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
   * Cold cache: a new {@link NamedAlleleMatcher} (and thus a new, empty {@code m_indexCache}) is created for every
   * sample, so the per-gene {@link HaplotypeCandidateIndex} is rebuilt from scratch every call.
   */
  @Benchmark
  public void coldCache(Blackhole bh) throws Exception {
    NamedAlleleMatcher matcher = new NamedAlleleMatcher(TestVcfBuilder.DEFAULT_TEST_ENV, m_definitionReader,
        /*findCombinations=*/false, /*topCandidateOnly=*/true, /*callCyp2d6=*/false);
    Result result = matcher.call(m_vcfFile, null, null);
    bh.consume(result);
  }


  /**
   * Warm cache: the same {@link NamedAlleleMatcher} instance is reused across samples, so the per-gene
   * {@link HaplotypeCandidateIndex} is built once and shared thereafter.
   */
  @Benchmark
  public void warmCache(Blackhole bh) throws Exception {
    Result result = m_warmMatcher.call(m_vcfFile, null, null);
    bh.consume(result);
  }
}
