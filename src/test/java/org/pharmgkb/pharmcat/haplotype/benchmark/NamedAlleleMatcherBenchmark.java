package org.pharmgkb.pharmcat.haplotype.benchmark;

import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.pharmgkb.pharmcat.Env;
import org.pharmgkb.pharmcat.ReportableException;
import org.pharmgkb.pharmcat.VcfFile;
import org.pharmgkb.pharmcat.definition.DefinitionReader;
import org.pharmgkb.pharmcat.haplotype.NamedAlleleMatcher;
import org.pharmgkb.pharmcat.util.DataManager;


/**
 * End-to-end benchmark harness for {@link NamedAlleleMatcher#call}.
 *
 * <p>Excluded from the default {@code test} task via the {@code benchmark} JUnit tag; run with {@code ./gradlew
 * benchmark}. Reports min / median / std-dev / cv% wall-clock time per scenario. Tune iterations with the system
 * properties {@code benchmark.warmup} (default 30) and {@code benchmark.iterations} (default 100).</p>
 *
 * <h2>Tier 1 noise hardening</h2>
 * <p>With {@code warmup=30} and {@code iterations=100}, hot methods reliably reach HotSpot C2 during warmup,
 * shrinking the end-to-end median noise floor to roughly 1–2%. The coefficient-of-variation (cv%) column makes
 * the remaining noise visible at a glance.</p>
 *
 * <p>To further reduce OS-scheduler jitter, pin the JVM to a single core:</p>
 * <pre>
 *   taskset -c 0 ./gradlew benchmark
 * </pre>
 */
@Tag("benchmark")
public class NamedAlleleMatcherBenchmark {
  private static final int DEFAULT_WARMUP = 30;
  private static final int DEFAULT_ITERATIONS = 100;


  @Test
  void runBenchmarks(@TempDir Path workDir) throws Exception {
    int warmup = intProp("benchmark.warmup", DEFAULT_WARMUP);
    int iterations = intProp("benchmark.iterations", DEFAULT_ITERATIONS);
    boolean timing = Boolean.parseBoolean(System.getProperty("benchmark.timing", "false"));

    Env env = new Env();
    PositionsIndex positionsIndex = PositionsIndex.getInstance();

    List<Scenario> scenarios = Scenarios.all();
    Result[] results = new Result[scenarios.size()];

    System.out.println();
    System.out.println("NamedAlleleMatcher benchmark: " + scenarios.size() + " scenarios, warmup=" + warmup +
        ", iterations=" + iterations);
    System.out.println();

    for (int i = 0; i < scenarios.size(); i += 1) {
      Scenario scenario = scenarios.get(i);
      results[i] = runScenario(scenario, workDir, env, positionsIndex, warmup, iterations);
    }

    printReport(results);

    if (timing) {
      System.out.println();
      System.out.println("=== NamedAlleleMatcher per-stage timing ===");
      for (Scenario scenario : scenarios) {
        System.out.println();
        System.out.println("--- timing: " + scenario.name() + " ---");
        runTimingPass(scenario, workDir, env, positionsIndex);
      }
    }
  }


  private Result runScenario(Scenario scenario, Path workDir, Env env, PositionsIndex positionsIndex,
      int warmup, int iterations) throws IOException, ReportableException {

    DefinitionReader definitionReader = new DefinitionReader(
        List.of(DataManager.getDefinitionFilePath(scenario.gene())),
        DataManager.DEFAULT_EXEMPTIONS_FILE);

    BenchmarkVcfBuilder builder = new BenchmarkVcfBuilder(scenario.gene(), positionsIndex, definitionReader);
    scenario.setup().accept(builder);

    Path vcfPath = workDir.resolve(sanitize(scenario.name()) + ".vcf");
    builder.write(vcfPath);

    NamedAlleleMatcher matcher = new NamedAlleleMatcher(env, definitionReader,
        scenario.findCombinations(), scenario.topCandidateOnly(), scenario.callCyp2d6());

    VcfFile sharedVcfFile = new VcfFile(vcfPath);
    // Warmup both patterns so JIT state is stable before either measurement pass.
    for (int w = 0; w < warmup; w += 1) {
      matcher.call(new VcfFile(vcfPath), null);
      matcher.call(sharedVcfFile, null);
    }

    // e2e: fresh VcfFile per iteration -> includes file I/O + VcfReader parse + matcher.
    long[] e2eNanos = new long[iterations];
    for (int m = 0; m < iterations; m += 1) {
      long t0 = System.nanoTime();
      matcher.call(new VcfFile(vcfPath), null);
      e2eNanos[m] = System.nanoTime() - t0;
    }

    // match: shared VcfFile -> VcfFile caches file bytes after first open, so file I/O is skipped.
    // VcfReader still parses each call, so this is not pure matcher cost, but the delta vs. e2e isolates
    // the file-open + read-all-bytes overhead.
    long[] matchNanos = new long[iterations];
    for (int m = 0; m < iterations; m += 1) {
      long t0 = System.nanoTime();
      matcher.call(sharedVcfFile, null);
      matchNanos[m] = System.nanoTime() - t0;
    }

    return new Result(scenario.name(), e2eNanos, matchNanos);
  }


  private void runTimingPass(Scenario scenario, Path workDir, Env env, PositionsIndex positionsIndex)
      throws IOException, ReportableException {

    DefinitionReader definitionReader = new DefinitionReader(
        List.of(DataManager.getDefinitionFilePath(scenario.gene())),
        DataManager.DEFAULT_EXEMPTIONS_FILE);

    BenchmarkVcfBuilder builder = new BenchmarkVcfBuilder(scenario.gene(), positionsIndex, definitionReader);
    scenario.setup().accept(builder);

    Path vcfPath = workDir.resolve(sanitize(scenario.name()) + "_timing.vcf");
    builder.write(vcfPath);

    NamedAlleleMatcher matcher = new NamedAlleleMatcher(env, definitionReader,
        scenario.findCombinations(), scenario.topCandidateOnly(), scenario.callCyp2d6()).timing();

    matcher.call(new VcfFile(vcfPath), null);
  }


  private static void printReport(Result[] results) {
    int nameWidth = Arrays.stream(results).mapToInt(r -> r.name.length()).max().orElse(20);
    nameWidth = Math.max(nameWidth, 20);

    String header = String.format(Locale.ROOT,
        "%-" + nameWidth + "s  %10s  %10s  %10s  %7s  %10s  %10s  %10s  %7s  %8s",
        "scenario", "e2e min", "e2e med", "e2e sd", "e2e cv%",
        "match min", "match med", "match sd", "match cv%", "delta%");
    String sep = "-".repeat(header.length());
    System.out.println();
    System.out.println("=== NamedAlleleMatcher benchmark results (ms) ===");
    System.out.println("sd = population std-dev (ms);  cv% = coefficient of variation = std-dev / mean");
    System.out.println("e2e:   fresh VcfFile per iteration (file I/O + VcfReader parse + matcher)");
    System.out.println("match: shared VcfFile (VcfReader parse + matcher, no file I/O)");
    System.out.println("delta%: fraction of e2e median attributable to file I/O");
    System.out.println(sep);
    System.out.println(header);
    System.out.println(sep);
    for (Result r : results) {
      double e2eMedian = r.medianMs(r.e2eNanos);
      double matchMedian = r.medianMs(r.matchNanos);
      double deltaPct = e2eMedian > 0 ? (e2eMedian - matchMedian) / e2eMedian * 100.0 : 0.0;
      System.out.printf(Locale.ROOT,
          "%-" + nameWidth + "s  %10.3f  %10.3f  %10.3f  %6.1f%%  %10.3f  %10.3f  %10.3f  %6.1f%%  %7.1f%%%n",
          r.name,
          r.minMs(r.e2eNanos), e2eMedian, r.stdDevMs(r.e2eNanos), r.cvPct(r.e2eNanos),
          r.minMs(r.matchNanos), matchMedian, r.stdDevMs(r.matchNanos), r.cvPct(r.matchNanos),
          deltaPct);
    }
    System.out.println(sep);
  }


  private static int intProp(String key, int fallback) {
    String v = System.getProperty(key);
    if (v == null || v.isEmpty()) {
      return fallback;
    }
    return Integer.parseInt(v);
  }

  private static String sanitize(String s) {
    return s.replaceAll("[^A-Za-z0-9_.-]+", "_");
  }


  private record Result(String name, long[] e2eNanos, long[] matchNanos) {
    double minMs(long[] nanos) {
      return Arrays.stream(nanos).min().orElse(0) / 1_000_000.0;
    }
    double medianMs(long[] nanos) {
      long[] sorted = nanos.clone();
      Arrays.sort(sorted);
      int n = sorted.length;
      long v = (n % 2 == 1) ? sorted[n / 2] : (sorted[n / 2 - 1] + sorted[n / 2]) / 2;
      return v / 1_000_000.0;
    }
    double meanMs(long[] nanos) {
      return Arrays.stream(nanos).average().orElse(0) / 1_000_000.0;
    }
    double stdDevMs(long[] nanos) {
      double mean = Arrays.stream(nanos).average().orElse(0);
      double variance = Arrays.stream(nanos).mapToDouble(x -> (x - mean) * (x - mean)).average().orElse(0);
      return Math.sqrt(variance) / 1_000_000.0;
    }
    double cvPct(long[] nanos) {
      double mean = meanMs(nanos);
      return mean > 0 ? stdDevMs(nanos) / mean * 100.0 : 0.0;
    }
  }
}
