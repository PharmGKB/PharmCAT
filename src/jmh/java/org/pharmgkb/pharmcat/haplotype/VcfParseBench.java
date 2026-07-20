package org.pharmgkb.pharmcat.haplotype;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.io.InputStreamReader;
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
import org.pharmgkb.parser.vcf.VcfParser;
import org.pharmgkb.parser.vcf.model.VcfMetadata;
import org.pharmgkb.parser.vcf.model.VcfPosition;
import org.pharmgkb.parser.vcf.model.VcfSample;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * JMH microbenchmark sizing the per-sample VCF read/parse cost that {@code NamedAlleleMatcher.call()} pays on every
 * sample, using the real all-gene single-sample {@code reference.vcf} (1226 data lines).
 *
 * <p>Run with: {@code ./gradlew jmh}</p>
 *
 * <ul>
 *   <li>{@code parseFull} — the full per-sample path: {@link VcfFile#getReader}.getAlleleMap(), i.e. vcf-parser
 *       tokenization + {@code VcfReader.parseLine} interpretation. This is what {@code call()} actually pays.</li>
 *   <li>{@code tokenizeOnly} — the same file run through {@code VcfParser} with a no-op line parser, isolating the
 *       vcf-parser library's tokenization + model construction (VcfPosition + every VcfSample).</li>
 * </ul>
 *
 * <p>The difference {@code parseFull - tokenizeOnly} approximates PharmCAT's own per-line work, which tells us whether
 * to optimize the vcf-parser library, PharmCAT's {@code VcfReader}, or both. {@code VcfFile} reads the file into memory,
 * so both benchmarks measure parse CPU without per-invocation disk I/O.</p>
 */
@State(Scope.Benchmark)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 5, time = 2)
@Measurement(iterations = 5, time = 2)
@Fork(2)
public class VcfParseBench {

  private Path m_tempDir;
  private Path m_vcfPath;
  private VcfFile m_vcfFile;
  private DefinitionReader m_definitionReader;
  private byte[] m_vcfBytes;


  @Setup(Level.Trial)
  public void setup() throws Exception {
    m_tempDir = Files.createTempDirectory("jmh-vcfparse-");
    m_vcfPath = m_tempDir.resolve("reference.vcf");
    try (InputStream in = VcfParseBench.class.getResourceAsStream("/org/pharmgkb/pharmcat/reference.vcf")) {
      if (in == null) {
        throw new IllegalStateException("Could not find reference.vcf on the classpath");
      }
      Files.copy(in, m_vcfPath);
    }
    m_vcfBytes = Files.readAllBytes(m_vcfPath);

    // full definition set (all genes) so locationsOfInterest matches a real run
    m_definitionReader = new DefinitionReader(DataManager.DEFAULT_DEFINITION_DIR, null,
        DataManager.DEFAULT_EXEMPTIONS_FILE);
    m_vcfFile = new VcfFile(m_vcfPath);
    // prime the in-memory buffer and sample list
    m_vcfFile.getReader(m_definitionReader, null, false).getAlleleMap();
  }


  @TearDown(Level.Trial)
  public void tearDown() throws Exception {
    Files.deleteIfExists(m_vcfPath);
    Files.deleteIfExists(m_tempDir);
  }


  @Benchmark
  public void parseFull(Blackhole bh) throws Exception {
    bh.consume(m_vcfFile.getReader(m_definitionReader, null, false).getAlleleMap());
  }


  @Benchmark
  public void tokenizeOnly(Blackhole bh) throws Exception {
    CountingLineParser counter = new CountingLineParser();
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(m_vcfBytes)));
         VcfParser parser = new VcfParser.Builder().fromReader(reader).parseWith(counter).build()) {
      parser.parseMetadata();
      parser.parse();
    }
    bh.consume(counter.count);
  }


  /** Consumes every parsed line without doing any PharmCAT-specific interpretation. */
  private static final class CountingLineParser implements org.pharmgkb.parser.vcf.VcfLineParser {
    private int count;

    @Override
    public void parseLine(VcfMetadata metadata, VcfPosition position, List<VcfSample> sampleData) {
      count += 1;
    }
  }
}
