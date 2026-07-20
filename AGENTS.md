# AGENTS.md

This repository contains PharmCAT, a Java-based pharmacogenomics tool with a Python preprocessor and a small Node/Yarn
release/docs toolchain.

## Scope and intent

- Keep changes narrowly scoped to the user’s request.
- Prefer minimal edits that preserve existing behavior outside the requested area.
- Do not update generated data, examples, release metadata, or packaged artifacts unless the task explicitly requires it.

## Repository layout

- `src/`: main Java application code and resources.
- `preprocessor/`: Python-based preprocessing scripts and support files.
- `docs/`: documentation site content.
- `bin/`: packaged entrypoint scripts.
- `dockstore/`: pipeline packaging assets.
- `.github/workflows/`: CI and release workflows.

## Toolchain

- Java 17 is the baseline runtime and compile target.
- Gradle is the primary build/test entrypoint via `gradlew.bat` on Windows and `./gradlew` on Unix.
- Node `> 24` with Yarn 4 is used for release automation and docs deployment tasks, not the main application build.
- Python is used under `preprocessor/`, but Python package changes should stay isolated to that area unless explicitly
  requested.

## Common commands

From the repository root:

- Build/test Java code:
  - `gradlew.bat test`
  - `gradlew.bat shadowJar`
- Test Python preprocessor:
  - `make test-preprocessor`
- Run Java tests with timing and a pass/fail summary (Unix helper; wraps `./gradlew cleanTest test`):
  - `bin/run-tests.sh` — full test suite
  - `bin/run-tests.sh 'org.pharmgkb.pharmcat.DpydTest' 'org.pharmgkb.pharmcat.haplotype.*'` — one or more Gradle
    `--tests` filter patterns
  - Prints wall-clock duration, the Gradle exit code, and `tests=… failures=… errors=… skipped=…` (plus any failing
    test names); full output is written to `build/run-tests.log`.

Prefer the smallest command that validates the change.

## Files that are commonly generated or release-managed

Be deliberate when editing these. Only change them when the task explicitly calls for it.

- `CHANGELOG.md`
- `docs/_config.yml`
- `bin/setup.sh`
- `preprocessor/pcat/common.py`
- `dockstore/pipeline/README.md`
- `dockstore/pipeline/PharmCAT_Pipeline.wdl`

`package.json` drives semantic-release configuration for the files above.
 
Never modify these files: 
- top-level `pharmcat_positions.*` and `pharmcat_regions.bed`


## Editing guidance

- Follow `.editorconfig`: spaces by default, 2-space indentation, 4 spaces for Python, tabs in `Makefile`.
- Preserve existing naming and structure conventions in touched files.
- Do not introduce unrelated formatting churn.
- For Java changes, keep compatibility with the current Gradle/JDK 17 setup unless the task explicitly changes toolchain
  requirements.

## Validation guidance

- For Java logic changes, prefer targeted Gradle tests first (e.g. `bin/run-tests.sh '<test pattern>'...`, which
  reports a pass/fail summary).
- For preprocessor-only changes, validate the smallest relevant command or test path in `preprocessor/`.
- For docs-only changes, avoid runtime-heavy validation unless the task depends on generated docs output.
- If a task touches release or packaging behavior, note clearly which outputs are expected to change.

## Cautions

- This is a biomedical/pharmacogenomics project. Treat behavioral changes as high sensitivity; avoid silent logic
  changes without corresponding validation.
- Some workflows update generated examples, caches, and reference assets. Do not run them casually in the course of
  unrelated work.
- The `Makefile` includes commands that clean and regenerate test data/results. Use those only when the task requires
  that workflow.
