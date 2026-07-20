# Benchmarking `NamedAlleleMatcher`

Two complementary tools live in the repo. They answer different questions.

- **Tier 1 — end-to-end harness** (JUnit): measures wall-clock time of `NamedAlleleMatcher.call()`
  across representative gene/mode scenarios. Use when the question is user-facing performance.
- **Tier 2 — JMH microbenchmarks**: measures individual hot loops in nanoseconds with proper
  warmup, forking, and statistical output. Use when the question is whether a small
  inner-loop change actually helps.

Neither replaces the other. Tier 1 tells you what users see; Tier 2 tells you whether a
targeted optimization is real.


## Tier 1: end-to-end harness

**Location**: `src/test/java/org/pharmgkb/pharmcat/haplotype/benchmark/`

- `NamedAlleleMatcherBenchmark.java` — the JUnit test, tagged `benchmark` so it's excluded
  from the default `test` task.
- `Scenarios.java` — catalog of 23 scenarios covering every `NamedAlleleMatcher` code path
  (standard/all-ref/single-het/multi-het, combinations mode, lowest-function-gene stages,
  DPYD HapB3 variants, missing positions).
- `BenchmarkVcfBuilder.java` + `PositionsIndex.java` — generate scenario VCFs from
  `pharmcat_positions.vcf` so scenarios stay in sync with the canonical position list.

### Running

```
./gradlew benchmark
```

Pinned to one core (reduces OS-scheduler jitter, recommended for comparative runs):

```
taskset -c 0 ./gradlew benchmark
```

Tune via system properties (defaults in parens):

- `benchmark.warmup` (30) — JVM warmup iterations per scenario before measuring.
- `benchmark.iterations` (100) — measured iterations per scenario.
- `benchmark.timing` (false) — after the summary, run each scenario once more with the
  built-in `MatcherTimings` instrumentation enabled so per-stage timings print alongside
  the wall-clock table.

Example: `./gradlew benchmark -Dbenchmark.iterations=250 -Dbenchmark.timing=true`.

### Report format

Two measurement passes per scenario:

- **e2e** — fresh `VcfFile` per iteration. Includes file I/O + `VcfReader` parse + matcher.
  Matches what users experience.
- **match** — shared `VcfFile` (bytes cached after first read). Skips file I/O, so the
  delta between e2e and match isolates parse+I/O overhead.

Columns per pass: `min`, `median`, `std-dev`, `cv%` (coefficient of variation =
std-dev / mean). `cv%` is the "noise floor" indicator: if it's larger than the delta
you're chasing, the delta isn't real.

### Limitations

1. **Noise floor is workload-dependent.** Large workloads (RYR1 unphased hets, CYP2D6
   combinations) settle around `cv%` 3–10%. Small workloads (< 3 ms) can be 30–70% even
   at 30/100. Trust `cv%` — don't chase deltas below it.

2. **Cross-scenario JIT contamination.** All 23 scenarios share one JVM. Later scenarios
   inherit JIT state from earlier ones. Reordering `Scenarios.all()` can shift numbers.
   The current design accepts this because it keeps the harness cheap and the shared code
   paths (which are what we care about) all reach C2 well before measurement starts.

3. **VcfReader parse is baked in.** Even the "match" column includes `VcfReader.parse`
   because `matcher.call()` re-parses on every call. If a matcher-internal change is
   small relative to parse cost, it will get diluted here.

4. **Sub-microsecond changes are invisible.** A change that saves 10 ns/call on an inner
   loop called 10k times per scenario changes wall clock by 100 μs. That's below the
   noise floor for every scenario. Reach for Tier 2 instead.

5. **Single machine, single run per invocation.** No independent JVM per scenario, no
   automatic replication. For before/after comparisons, run baseline and change back-to-back
   in the same shell session with the same load conditions.


## Tier 2: JMH microbenchmarks

**Location**: `src/jmh/java/org/pharmgkb/pharmcat/haplotype/`

Handled by the `me.champeau.jmh` Gradle plugin. Configuration in `build.gradle`:

```groovy
jmh {
  warmupIterations = 5
  iterations = 5
  fork = 2
  resultFormat = 'TEXT'
  includeTests = true   // JMH source can see test helpers
}
```

### Running

```
./gradlew jmh
```

Each `@Benchmark` runs 2 forks × 5 warmup × 5 measurement = 10 samples with C2-warm JIT.
Output includes mean and standard error; non-overlapping CIs are the go/no-go signal.

### Existing benchmarks

- `CoreMatchingPathBench.java` — the per-sample core matching path (`MatchData` construction,
  `marshallHaplotypes`, `generateSamplePermutations`, shared `HaplotypeCandidateIndex` attach, and
  `comparePermutations`) for CYP2C19 and RYR1, in phased (≤2 permutations) and unphased (2^n permutations) regimes.
- `ColdVsWarmIndexCacheBench.java` — cold `HaplotypeCandidateIndex` build vs. the warm shared-cache lookup.

### Adding a new benchmark

1. Put the class in `src/jmh/java/org/pharmgkb/pharmcat/haplotype/` if you need access to
   package-private members of `MatchData`, `DiplotypeMatcher`, `SamplePermutation`, etc.
   Otherwise use a subpackage.
2. Follow the pattern in `CoreMatchingPathBench`:
   - `@State(Scope.Benchmark)` with `@Setup(Level.Trial)` that builds a realistic
     `MatchData` via the same call sequence `NamedAlleleMatcher.initializeCallData` uses.
   - One `@Benchmark` per case you want reported side by side — per scenario (as
     `CoreMatchingPathBench` does for gene × phased/unphased) or, when comparing candidate
     implementations, one per variant so JMH measures them in the same run.
   - Return a value from each `@Benchmark` so JMH's DCE guard doesn't strip work.
3. Reuse `Scenarios` / `BenchmarkVcfBuilder` / `PositionsIndex` where possible — they're
   on the JMH classpath via `includeTests = true`.

### Limitations

1. **You have to reproduce production state.** JMH benchmarks build objects directly
   rather than going through `NamedAlleleMatcher.call()`. Get the setup wrong and you're
   measuring the wrong thing. Study the real init path before mocking it.

2. **Microbenchmarks lie by omission.** A 2.5× win in an inner loop that runs for 100 μs
   per matcher call saves at most 60 μs total. If that's not visible in Tier 1, it may
   still be worth doing (correctness, readability) but don't oversell the user-facing
   impact.

3. **JMH itself has fixed overhead.** Expect a few ns of framework overhead per
   `@Benchmark` call. Fine for anything ≥ 5 ns/op; misleading below that.


## Decision guide

| Question | Tool |
|---|---|
| "Did my change make PharmCAT faster for users?" | Tier 1, look at `match` column |
| "Is `VcfReader` parse a meaningful fraction of matcher cost?" | Tier 1, `e2e - match` delta |
| "Where does the time go inside `matcher.call()` for scenario X?" | Tier 1 with `-Dbenchmark.timing=true` |
| "Did this change to hot loop Y actually help?" | Tier 2, JMH |
| "Which of these two candidate implementations is faster?" | Tier 2, JMH — one `@Benchmark` per candidate |
| "Is a sub-1% end-to-end delta real or noise?" | Tier 2 on the specific loop that changed |

