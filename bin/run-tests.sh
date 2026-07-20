#!/usr/bin/env bash
#
# Dev helper: run Gradle tests with timing, exit-code capture, and a pass/fail summary.
# Forces a clean re-run (cleanTest) so the reported counts reflect only this invocation.
#
# Usage:
#   eng/run-tests.sh                          # full test suite
#   eng/run-tests.sh <pattern>...             # one or more Gradle --tests filter patterns, e.g.
#   eng/run-tests.sh 'org.pharmgkb.pharmcat.DpydTest' 'org.pharmgkb.pharmcat.haplotype.*'
#
set -uo pipefail

cd "$(dirname "$0")/.." || exit 2
mkdir -p build
log="build/run-tests.log"

args=()
for pat in "$@"; do
  args+=(--tests "$pat")
done

start=$(date +%s)
./gradlew cleanTest test "${args[@]}" --console=plain > "$log" 2>&1
exit_code=$?
end=$(date +%s)
dur=$(( end - start ))

echo "duration: ${dur}s ($(( dur / 60 ))m $(( dur % 60 ))s), gradle exit: ${exit_code}"
grep -E "BUILD SUCCESSFUL|BUILD FAILED" "$log" | tail -1 || true

python3 - <<'PY'
import glob, xml.etree.ElementTree as ET
t = f = e = s = 0
fails = []
for path in glob.glob("build/test-results/test/*.xml"):
    r = ET.parse(path).getroot()
    t += int(r.get("tests", 0)); f += int(r.get("failures", 0))
    e += int(r.get("errors", 0)); s += int(r.get("skipped", 0))
    for tc in r.findall("testcase"):
        if tc.findall("failure") or tc.findall("error"):
            fails.append(r.get("name", "").split(".")[-1] + "." + tc.get("name", ""))
print(f"tests={t} failures={f} errors={e} skipped={s}")
if fails:
    print("FAILING: " + ", ".join(fails[:30]))
PY

echo "full log: $log"
exit "$exit_code"
