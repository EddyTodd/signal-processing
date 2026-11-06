# Numerical validation and research-grade evidence

This library is actively evaluated; passing included tests does not prove global error bounds or superiority to external DSP libraries.

## Validation tiers

1. Independent textbook O(N²) binary64 DFT and linear-convolution oracles.
2. Parseval's theorem, circular-shift phase, inverse round-trip and DCT-II equivalence.
3. Existing algorithm and filter design regression suites.
4. Independent AddressSanitizer and UndefinedBehaviorSanitizer build where available.
5. External Local CI microbenchmarks with compiler, flags, CPU, ISA, samples, spread, repetitions and allocation regime documented.

## Crucial Release-build contract

The test binaries use C assert statements. A normal Release configuration defines NDEBUG, disabling those assertions. CMake now explicitly undefines NDEBUG on test binaries for every build type, and the numerical test deliberately fails compilation if this policy is violated.

## Limitations

Textbook binary64 oracles are algorithmically independent, but not exact arbitrary-precision reference calculations. Finite seeded tests are regression evidence, not formal proof. Performance data must be measured and attributed to a particular machine and build. Never present illustrative numbers as measured results.

## Measured baseline — 2026-11-06

A full Local CI campaign on source commit `dd299dacc2a4b4ef2c01953c6304da2c8a27de0c`
passed eight Release CTest executables and eight ASan/UBSan CTest executables.
A separate single-threaded binary64 campaign produced 118 algorithm/shape
cases with nine repetitions per case, 1,062 raw timing samples total, and
reported medians, empirical p10–p90, and median absolute deviations.
The machine was a 13th Gen Intel Core i7-13700KF running WSL with GCC 13.3.0.
This is a same-machine timing baseline, not a claim of optimized code,
universal results, statistically representative hardware, or peer review.
Exact GitHub report: https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/dd299dacc2a4-bdccb989-1d39f279acb77f33

## Independent extended-precision tests

The additional `numerical_precision` suite uses a separate O(N²) direct DFT
accumulated in `long double` and tests both `float` and `double` kernels
for prime/composite and radix-two lengths, as well as impulse, constant,
alternating, coherent-tone, random, and mixed-dynamic-range inputs.
It includes multiple FFT families and inverse DFT reconstruction checks.
This suite is **pending actual CI** for the new source revision. Long-double
arithmetic is a higher-precision independent oracle on the tested platform,
not an arbitrary-precision guarantee for every architecture.

## Verified extended-precision baseline (2026-11-06)

Source `5e8707acd4b836828418f4177a4c39aa83b0dbc4`, Local CI run `5e8707acd4b8-358bd155-9b69a8ab6020425d`: all nine Release and all nine
ASan/UBSan CTest executables passed, including `numerical_precision` with
independent long-double O(N²) forward/inverse DFT comparison across float32
and float64 algorithms, structured signals and varied lengths. The 118-case
single-machine binary64 benchmark campaign also passed with nine timing
samples per case. This validates tested finite cases only; the 24-size,
additional-seed suite added after this run requires separate CI evidence.
