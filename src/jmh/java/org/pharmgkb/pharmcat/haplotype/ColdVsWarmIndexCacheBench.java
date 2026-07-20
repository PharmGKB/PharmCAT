package org.pharmgkb.pharmcat.haplotype;

import java.util.Arrays;
import java.util.List;
import java.util.SortedSet;
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
import org.openjdk.jmh.annotations.Warmup;
import org.openjdk.jmh.infra.Blackhole;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.HaplotypeCandidateIndex;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * JMH microbenchmark quantifying the payoff of issue #1: caching the per-gene {@link HaplotypeCandidateIndex} on
 * {@link DefinitionFile} so it is shared across all samples instead of rebuilt for every sample.
 *
 * <p>Run with: {@code ./gradlew jmh}</p>
 *
 * <p>{@code coldBuild} constructs a fresh {@link HaplotypeCandidateIndex} on every invocation, mirroring the
 * pre-#1 behavior where each sample rebuilt the index from scratch (the matcher's index cache lived on the
 * per-sample {@code NamedAlleleMatcher}). {@code warmLookup} fetches the index from the {@link DefinitionFile}
 * cache introduced by #1, mirroring the shared cross-sample (and cross-thread {@code BatchPharmCAT}) behavior. The
 * delta is the per-sample, per-gene cost that #1 eliminates.</p>
 *
 * <p>This is deliberately {@code Env}-free: the index depends only on the gene's callable {@link NamedAllele} set,
 * its sorted positions, and the {@code assumeReference} flag, none of which involve phenotype/drug data.</p>
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 2)
@Measurement(iterations = 5, time = 2)
@Fork(2)
public class ColdVsWarmIndexCacheBench {

  /** Gene whose definition drives the index (RYR1: many haplotypes/positions, standard reference-defaulted path). */
  private static final String GENE = "RYR1";

  private DefinitionFile m_definitionFile;
  private SortedSet<NamedAllele> m_haplotypes;
  private VariantLocus[] m_permutationPositions;


  @Setup(Level.Trial)
  public void setup() throws Exception {
    DefinitionReader definitionReader = new DefinitionReader(
        List.of(DataManager.getDefinitionFilePath(GENE)),
        DataManager.DEFAULT_EXEMPTIONS_FILE);
    m_definitionFile = definitionReader.getDefinitionFile(GENE);
    // Inputs to a cold build, matching DefinitionFile.buildCandidateIndex(assumeReference=true, findCombinations=false).
    m_haplotypes = m_definitionFile.getNamedAlleles();
    m_permutationPositions = Arrays.stream(m_definitionFile.getVariants()).sorted().toArray(VariantLocus[]::new);
    // Prime the warm cache once so warmLookup measures a pure cache hit.
    m_definitionFile.getCandidateIndex(true, false);
  }


  /**
   * Cold: rebuild the per-gene index from scratch every call (pre-#1 per-sample cost).
   */
  @Benchmark
  public void coldBuild(Blackhole bh) {
    bh.consume(new HaplotypeCandidateIndex(m_haplotypes, m_permutationPositions, true));
  }


  /**
   * Warm: fetch the index from the {@link DefinitionFile} cache introduced by #1 (shared across samples).
   */
  @Benchmark
  public void warmLookup(Blackhole bh) {
    bh.consume(m_definitionFile.getCandidateIndex(true, false));
  }
}
