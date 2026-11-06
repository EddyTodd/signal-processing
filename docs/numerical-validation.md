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

## Historical failed float32 size campaign (before the correction)

The larger 24-size/multi-seed suite failed 1/9 Release CTest executables on
`a07e7ec3aba685b76e2e6fb71cb561641ced8e7f`, while eight older suites
passed. An assertion reported a float32 error greater than the original
`0.0008 * (1 + |reference|)` bound, but omitted algorithm and specimen identity.
The follow-up adds explicit method/size/pattern/seed/phase/index and observed
error/bound diagnostics, **without changing the numerical acceptance tolerance**.
A subsequent diagnostic identified the direct-DFT DC cancellation defect; the verified correction is recorded below.

## Diagnosed float32 direct DFT cancellation (2026-11-06)

The exact CI diagnostic run `586e4b52f74d-348e2e17-e28b3d2a64762719`
identified the first expanded-suite failure as `float32 / direct / N=256 /
constant DC / k=224`, with observed absolute error ≈0.001234 and an
unchanged allowed bound of 0.0008. The independent long-double reference is
approximately zero in that bin. The direct O(N²) implementation used float32
angles, trigonometric roots and summation; this is a cancellation-sensitive
path. The proposed fix computes *only the float32 direct DFT internal phase
and accumulation* in binary64 and rounds to float32 at the output, while
retaining the public API and the binary64 path. It has **not yet passed local
CI** at the time of this commit. The subsequent pinned campaign below
verified finite-case accuracy and binary64 performance. The finite
error bound was not loosened.

## Verified expanded numerical campaign and direct DFT correction

Local CI report (2026-11-06):
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/19047b0782bf-32ac4220-d21dde9956499eed

Source `19047b0782bf03e3fc4e9716d73f25781c4ef590` with coordinator
`32ac4220897e3d41b16324f29eb682da0b65cc02` passed **9/9 Release**
and **9/9 ASan/UBSan** CTest executables (assertions enabled), covering the
expanded 24-length, multi-seed, float32/float64, forward/inverse independent
long-double O(N²) transform comparisons. The observed float32 direct DFT
constant-sequence cancellation failure is resolved for this finite suite
without loosening the prior 0.0008 × (1 + magnitude) threshold.

The associated i7-13700KF / GCC 13.3.0 binary64 benchmark artifact contains
**118** measured method/shape cases × **9** independent timing samples
(1,062 total), with median, p10–p90 empirical spread and MAD.
The website feature branch `b5a64e58ebe4e44d32507f720c58ff0ded6440a9`
now cites this exact published GitHub artifact, not synthetic values.

**Pending next Local CI:** A stricter analytic DC cancellation test now
exercises lengths 128, 256 and 257 in both transform directions. Separately,
coordinator branch `ai/portable-photocraft-coordinator` now builds both
binary32 and binary64 method templates for independent, matched-regime
microbenchmarks; fresh float32 samples and this new test must be reviewed
before reporting either verified. Frequency scaling/thermal state remain
uncontrolled, and no arbitrary-precision proof, cross-machine/ISA/compiler
or third-party FFTW/SciPy comparison has been completed.

## Matched binary64 direct-DFT timing check (two pinned campaigns)

Both artifacts used the *same* benchmark source hash, GCC 13.3.0,
i7-13700KF, `-O3 -DNDEBUG -std=c++23 -march=native`,
single-thread end-to-end allocating runs, and nine raw samples/case.

| N | Old median ns | New median ns | Old p10–p90 ns | New p10–p90 ns |
|---:|---:|---:|---:|---:|
| 16 | 2119.88 | 2103.40 | 2108.45–2182.51 | 2086.51–2227.86 |
| 32 | 9183.15 | 8921.06 | 9068.11–9432.31 | 8866.87–9033.13 |
| 64 | 38974.37 | 37662.78 | 37776.94–40330.42 | 37212.21–39006.37 |
| 128 | 182294.17 | 178681.20 | 180884.24–188237.14 | 176308.78–181511.51 |
| 256 | 690627.06 | 688272.75 | 684850.06–739604.93 | 685896.44–692804.73 |

Original source: `5e8707acd4b836828418f4177a4c39aa83b0dbc4`,
report `5e8707acd4b8-358bd155-9b69a8ab6020425d`.
New source: `19047b0782bf03e3fc4e9716d73f25781c4ef590`,
report `19047b0782bf-32ac4220-d21dde9956499eed`.
These are two distinct laptop runs; p10/p90 is within-run observed spread,
**not a confidence interval**. Differences do not establish a performance
improvement or regression, and these numbers describe the **binary64 path**,
not the newly corrected float32 path.
