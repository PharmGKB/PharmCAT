package org.pharmgkb.pharmcat.haplotype;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.SortedMap;
import java.util.SortedSet;
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
import org.pharmgkb.pharmcat.definition.model.DefinitionExemption;
import org.pharmgkb.pharmcat.definition.model.DefinitionFile;
import org.pharmgkb.pharmcat.definition.model.NamedAllele;
import org.pharmgkb.pharmcat.definition.model.VariantLocus;
import org.pharmgkb.pharmcat.haplotype.benchmark.BenchmarkVcfBuilder;
import org.pharmgkb.pharmcat.haplotype.benchmark.PositionsIndex;
import org.pharmgkb.pharmcat.haplotype.model.HaplotypeMatch;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * JMH microbenchmark measuring the per-sample core haplotype-matching path that {@code NamedAlleleMatcher} runs for a
 * standard gene: {@link MatchData} construction, {@code marshallHaplotypes}, {@code generateSamplePermutations}, the
 * shared {@link org.pharmgkb.pharmcat.definition.model.HaplotypeCandidateIndex} attach (issue #1), and
 * {@code comparePermutations}.
 *
 * <p>Run with: {@code ./gradlew jmh}</p>
 *
 * <p>This intentionally stops short of {@link DiplotypeMatcher} pairing/scoring, which requires an {@code Env} (and
 * therefore phenotype/drug data that is unrelated to matching and unavailable in the JMH runtime) and is unaffected by
 * #1. The sample VCF sets the first few positions heterozygous/unphased to force multiple permutations rather than a
 * trivial all-reference match. Each iteration re-runs the full per-sample pipeline, so the shared index attach reflects
 * the post-#1 cache-hit path.</p>
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 2)
@Measurement(iterations = 5, time = 2)
@Fork(2)
public class CoreMatchingPathBench {

  private GeneCase m_cyp2c19Unphased;
  private GeneCase m_cyp2c19Phased;
  private GeneCase m_ryr1Unphased;
  private GeneCase m_ryr1Phased;
  private Path m_tempDir;


  @Setup(Level.Trial)
  public void setup() throws Exception {
    m_tempDir = Files.createTempDirectory("jmh-corematch-");
    PositionsIndex positionsIndex = PositionsIndex.getInstance();
    m_cyp2c19Unphased = GeneCase.build("CYP2C19", positionsIndex, m_tempDir, false);
    m_cyp2c19Phased = GeneCase.build("CYP2C19", positionsIndex, m_tempDir, true);
    m_ryr1Unphased = GeneCase.build("RYR1", positionsIndex, m_tempDir, false);
    m_ryr1Phased = GeneCase.build("RYR1", positionsIndex, m_tempDir, true);
  }


  @TearDown(Level.Trial)
  public void tearDown() throws Exception {
    if (m_tempDir != null) {
      Files.deleteIfExists(m_cyp2c19Unphased.vcfPath);
      Files.deleteIfExists(m_cyp2c19Phased.vcfPath);
      Files.deleteIfExists(m_ryr1Unphased.vcfPath);
      Files.deleteIfExists(m_ryr1Phased.vcfPath);
      Files.deleteIfExists(m_tempDir);
    }
  }


  @Benchmark
  public void cyp2c19Unphased(Blackhole bh) {
    bh.consume(m_cyp2c19Unphased.run());
  }


  @Benchmark
  public void cyp2c19Phased(Blackhole bh) {
    bh.consume(m_cyp2c19Phased.run());
  }


  @Benchmark
  public void ryr1Unphased(Blackhole bh) {
    bh.consume(m_ryr1Unphased.run());
  }


  @Benchmark
  public void ryr1Phased(Blackhole bh) {
    bh.consume(m_ryr1Phased.run());
  }


  /** Precomputed per-gene inputs plus the run() that exercises the core matching path each invocation. */
  private static final class GeneCase {
    private final String gene;
    private final SortedMap<String, SampleAllele> alleleMap;
    private final VariantLocus[] allPositions;
    private final SortedSet<NamedAllele> haplotypes;
    private final @Nullable DefinitionExemption exemption;
    private final @Nullable SortedSet<VariantLocus> extraPositions;
    private final DefinitionFile definitionFile;
    private final Path vcfPath;

    private GeneCase(String gene, SortedMap<String, SampleAllele> alleleMap, VariantLocus[] allPositions,
        SortedSet<NamedAllele> haplotypes, @Nullable DefinitionExemption exemption, DefinitionFile definitionFile,
        Path vcfPath) {
      this.gene = gene;
      this.alleleMap = alleleMap;
      this.allPositions = allPositions;
      this.haplotypes = haplotypes;
      this.exemption = exemption;
      this.extraPositions = exemption == null ? null : exemption.getExtraPositions();
      this.definitionFile = definitionFile;
      this.vcfPath = vcfPath;
    }

    static GeneCase build(String gene, PositionsIndex positionsIndex, Path tempDir, boolean phased) throws Exception {
      DefinitionReader definitionReader = new DefinitionReader(
          List.of(DataManager.getDefinitionFilePath(gene)),
          DataManager.DEFAULT_EXEMPTIONS_FILE);
      BenchmarkVcfBuilder builder = new BenchmarkVcfBuilder(gene, positionsIndex, definitionReader);
      int n = Math.min(6, builder.size());
      // Phased het (0|1) yields 2 permutations (effectively phased); unphased het (0/1) yields 2^n permutations.
      String gt = phased ? "0|1" : "0/1";
      for (int i = 0; i < n; i++) {
        builder.set(i, gt);
      }
      Path vcfPath = tempDir.resolve(gene.toLowerCase() + (phased ? "_phased" : "_unphased") + ".vcf");
      builder.write(vcfPath);

      VcfFile vcfFile = new VcfFile(vcfPath);
      SortedMap<String, SampleAllele> alleleMap =
          vcfFile.getReader(definitionReader, null, false, false).getAlleleMap();
      return new GeneCase(gene, alleleMap, definitionReader.getPositions(gene),
          definitionReader.getHaplotypes(gene), definitionReader.getExemption(gene),
          definitionReader.getDefinitionFile(gene), vcfPath);
    }

    SortedSet<HaplotypeMatch> run() {
      MatchData data = new MatchData("sample", gene, alleleMap, allPositions, extraPositions, exemption);
      data.marshallHaplotypes(gene, haplotypes, false, true);
      data.generateSamplePermutations();
      if (data.getNumSampleAlleles() > 0 && data.getMissingPositions().isEmpty()) {
        data.useSharedIndex(definitionFile.getCandidateIndex(true, false));
      }
      return data.comparePermutations();
    }
  }
}
