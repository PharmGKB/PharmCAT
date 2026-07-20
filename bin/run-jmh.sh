#!/usr/bin/env bash
#
# Dev helper: build the JMH benchmark jar and run org.openjdk.jmh.Main with full CLI control
# (forks, iterations, -prof profilers, etc.). The me.champeau.jmh Gradle task can't pass through
# arbitrary JMH args, so this builds the jar + a runtime classpath and invokes JMH directly.
#
# Full output is written to build/run-jmh.log; a filtered summary (the final results table, plus
# any `-prof stack` RUNNABLE hot-method section) is printed to stdout, so callers never need to
# pipe or redirect the output themselves.
#
# Usage:
#   bin/run-jmh.sh <benchmark-regex> [jmh args...]
#
# Examples:
#   bin/run-jmh.sh CoreMatchingPathBench
#   bin/run-jmh.sh 'VcfParseBench' -f 2 -wi 3 -i 5 -w 1s -r 2s -prof gc
#   bin/run-jmh.sh 'VcfParseBench.parseFull' -f 1 -r 3s -prof stack
#
# All arguments are passed straight through to org.openjdk.jmh.Main.
#
set -uo pipefail

cd "$(dirname "$0")/.." || exit 2

if [[ $# -eq 0 ]]; then
  echo "usage: bin/run-jmh.sh <benchmark-regex> [jmh args...]" >&2
  exit 2
fi

mkdir -p build
log="build/run-jmh.log"
initScript="$(mktemp)"
trap 'rm -f "$initScript"' EXIT

cat > "$initScript" <<'GRADLE'
gradle.projectsEvaluated {
  rootProject.tasks.register('printJmhCp') {
    doLast {
      def files = new LinkedHashSet()
      def cfg = rootProject.configurations.findByName('jmhRuntimeClasspath')
      if (cfg != null) { files.addAll(cfg.files) }
      ['jmh', 'test', 'runtimeClasspath'].each { name ->
        def c = rootProject.configurations.findByName(name)
        if (c != null && c.canBeResolved) { files.addAll(c.files) }
      }
      ['build/classes/java/jmh', 'build/classes/java/main', 'build/classes/java/test',
       'build/resources/main', 'build/resources/test'].each { files.add(rootProject.file(it)) }
      rootProject.file('build/jmh-classpath.txt').text = files.collect { it.absolutePath }.join(':')
    }
  }
}
GRADLE

# Build the benchmark jar and (re)generate the runtime classpath in one invocation.
# Build output (incl. the JMH bytecode generator's chatter) goes to the log, not stdout.
if ! ./gradlew -q -I "$initScript" jmhJar printJmhCp > "$log" 2>&1; then
  echo "error: gradle build failed" >&2
  tail -25 "$log" >&2
  exit 2
fi

jar="$(ls build/libs/*-jmh.jar 2>/dev/null | head -1)"
if [[ -z "$jar" || ! -f build/jmh-classpath.txt ]]; then
  echo "error: could not locate JMH jar or classpath" >&2
  exit 2
fi

start=$(date +%s)
java -cp "${jar}:$(cat build/jmh-classpath.txt)" org.openjdk.jmh.Main "$@" > "$log" 2>&1
code=$?
end=$(date +%s)
dur=$(( end - start ))

# `-prof stack` hot methods (RUNNABLE thread state), if present.
awk '/\[Thread state: RUNNABLE\]/{f=1; print; next}
     /\[Thread state: (TIMED_WAITING|WAITING|BLOCKED|NEW|TERMINATED)\]/{f=0}
     f' "$log"

# Final results table (JMH prints it after "# Run complete"; includes -prof gc secondary rows).
sed -n '/# Run complete/,$p' "$log" | awk '/^Benchmark[[:space:]]+Mode/{c=1} c'

if [[ $code -ne 0 ]]; then
  echo "--- last 25 lines of $log ---"
  tail -25 "$log"
fi

echo
echo "duration: ${dur}s ($(( dur / 60 ))m $(( dur % 60 ))s), jmh exit: ${code}"
echo "full log: $log"
exit "$code"
