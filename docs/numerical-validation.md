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

## Independent IIR oracle extension — pending fresh Local CI

A separate CTest target `signal_processing_iir_oracles` now checks binary32
and binary64 IIR realizations against a direct textbook difference-equation
reference accumulated in `long double`. Its known conjugate poles span
radii 0.05–0.985 and distinct angular positions. Deterministic impulse, step,
alternating and seeded random stimuli test eight realization paths, reset and
block/scalar equivalence, plus span-length exceptions. Exact first-order
impulse-response formulas for positive/negative real poles and non-unit
leading denominator coefficients are checked independently. A separate complex
rational response oracle checks DC, Nyquist, interior frequencies, and known
stability-domain boundaries. These tests are **not yet verified** by Local CI;
no formal stability or conditioning guarantee is implied for all parameters.

## Verified 236-case dual-precision local campaign (2026-11-06)

Local CI source/handoff `60f0923364188f9cfb728c40676c1adf6015e063`
and coordinator `cc8a52e49ba0c06578c2ec698c0bc605fa40a1c5`:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/60f092336418-cc8a52e4-85de15eedeb8b224

- **9/9 Release** tests and **9/9 AddressSanitizer/UndefinedBehaviorSanitizer**
  tests passed with assertions explicitly enabled, including the strict
  N=128/256/257 analytic float32/float64 cancellation regression.
- **236 distinct method/shape/precision cases:** 118 separately compiled
  binary32, 118 binary64. Nine repeated steady-clock timing samples per case
  (2,124 raw samples). The published artifact SHA-256 is
  `e5eee46e34bdbb085bef75c79c759c5c8ac62b8b3d7c6365bcf71e0911d5b5f9`.
- Intel i7-13700KF, GCC 13.3.0, `-O3 -DNDEBUG -std=c++23 -march=native`,
  one benchmark thread; availability of CPU ISA flags is recorded, but
  instruction selection was not independently inspected.
- Illustrative **measured rows** from the published artifact for N=256,
  end-to-end allocating DFT (ns per call, median with empirical p10–p90):

| Precision | Method | Median ns | p10–p90 ns |
|---|---|---:|---:|
| binary32 | direct DFT | 726838.44 | 713907.51–746311.83 |
| binary32 | radix-2 | 5195.98 | 5135.20–5407.87 |
| binary32 | Bluestein | 40107.52 | 39399.67–42167.84 |
| binary64 | direct DFT | 689812.00 | 686099.49–735952.29 |
| binary64 | radix-2 | 5157.75 | 5095.11–5199.23 |
| binary64 | Bluestein | 42334.65 | 41752.45–42979.69 |

Do **not** read this table as a universal benchmark or conclude a speedup
between precisions. Measured comparisons must match precision, shape and
execution regime. CPU frequency/thermal state, inter-run noise, alternate
compilers, platforms and external reference libraries remain uncontrolled.
The website `feature/signal-processing-lab-20261105` commit
`096e02aaa61bf345dc7fc7e50549ddb9afd4711c` displays the entire
verified artifact with a precision selector and within-precision ratios.
Its GitHub workflow `42169512038` passed Astro type checking, the
three measured-data tests, the site build and smoke checks after a prior
nullability failure was corrected.

**Next correctness gate:** the 10th new IIR oracle CTest executable has not yet
run on the laptop. The 9/9 results above validate only the preceding source
revision. Keep the source task ACTIVE pending new Release/ASan/UBSan evidence.

## IIR oracle compilation blocker (2026-11-06)

Pinned failed run:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/d066ee6d0865-0d022547-ec5d160c55a751b2

Source `d066ee6d0865b7043a0b214ad384caa25a236930` failed at
`release-numerical-tests` while compiling the newly introduced
`tests/iir_oracles.cpp`. GCC 13 reported an **ambiguous brace-list
constructor** for `BiquadTransposedDirectFormII<T>` in the invalid-domain
test, for both float32 and float64. This is a *test harness compile error*;
no CTest numerical checks, ASan/UBSan, or timing campaign executed in this
run. The fix explicitly constructs `BiquadCoefficients<T>` as the argument
without changing the library or the invalid-`a0` expectation. Added
specimen/method-specific recurrence and frequency-response failure messages,
retaining the original tolerances. These fixes remain pending fresh Local CI.
Previous passing source `60f0923` remains the current verified 9/9 and
236-case dual-precision evidence.

## Second independently pinned 10-suite pass, and sparse-hop STFT research step

The user's latest report for source `6e1862fba245051f75205401f9ae37981cd533e1`, coordinator
`bff04cc91ded737fa5cd7803ad344126f9d613dc`, is:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/6e1862fba245-bff04cc9-fbded4472b5b6c9c

All **10/10 Release** and **10/10 ASan/UBSan** tests passed with assertions on,
including the separately compiled float32/float64 independent IIR reference.
The 236-case FFT/convolution/DCT/FIR native binary32/binary64 campaign reports
nine raw samples per case (2,124 total), empirical medians, p10/p90, MAD,
GCC 13.3.0, i7-13700KF, `-march=native`, CPU capability flags, and pinned
benchmark source hash `85a3c51d593800ec2c7fc137aeee1bc2f86701c04cf999d66a6bf0fef73a4198`.
Published benchmark-artifact SHA-256:
`4e2318c191c1e6cb15e12612af2d5876d3eee4d954b05e3167173497416ed122`.

The earlier `acfb166f` run uses the identical *source*, benchmark code,
compiler flags and machine, but the timing medians differ for all 236 cases.
They are two noisy runs, not evidence of an algorithmic speedup or slowdown.
The website's 236-case precision selector/chart/table now points only to
the latest real pinned run, in feature commit
`a8c1459e82c7933a515e7e4d0fd76c6c6149e0ec`.
GitHub Actions run
https://github.com/EddyTodd/eddytodd-web/actions/runs/42171664762
passed Astro type-checking, all three raw-sample provenance tests,
the static build and smoke checks. Dynamic table method names retain
semantic row-header cells.

**Not yet verified:** The next source commit `5a48eed796e5fcc8c403342fe90c8fb76696826d`
bounds padded STFT frame counts for hop larger than frame size (avoids
fully empty synthetic tail frames) and introduces an independent long-double
STFT/ISTFT/PSD/Welch/CSD suite. A separate exhaustive frame-grid consistency
check compared the frame-count formula with independent enumeration over
866,880 small parameter combinations and found no mismatches; this is
*not* a Release CTest, sanitizer result or proof for all size_t inputs.
The new source requires 11/11 Release and 11/11 ASan/UBSan validation
before reporting the boundary fix as tested.

## Verified sparse-hop STFT and independent spectral oracles (2026-11-06)

Pinned source/handoff `33d867573a6be76e0cfc9b2296b82029d22ce865`, coordinator
`15329221c616b001e86f4fe08bbfa9e1989c2123`,
run `33d867573a6b-15329221-e4e57b64052b98af`:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/33d867573a6b-15329221-e4e57b64052b98af

The **11/11 Release** and **11/11 ASan/UBSan** suites passed
with explicit `-UNDEBUG` assertions, including the new
`signal_processing_spectral_oracles` executable covering independent
`long double` direct DFT, real/complex binary32/binary64 STFT and inverse
weighted overlap, spectral density, Welch, CSD and end padding/hop gaps.
The sparse-hop STFT synthetic-empty-tail correction is included in the
tested revision. Finite tests cannot prove all boundary inputs or bounds.

The pinned artifact contains **236 independently compiled**
method/shape/precision combinations (118 binary32 + 118 binary64), each
with nine repeated timed samples (**2,124 total**), measured on the
i7-13700KF with GCC 13.3.0 and `-O3 -march=native`. The reported CPU
ISA flags indicate availability, not verified instructions used.
Empirical p10/p90/MAD describe observed variability, not uncertainty
from repeated independent machines. Artifact SHA-256:
`4969f0d6fab60db4b99ce20881178d2e92236afcb66cac323423e971f192ba5a`.

The public website feature commit
`1609463077c842de0869e2ad29373ef1067f51d9`
references this exact measured GitHub artifact, filters precision,
does not mix unmatched speedup ratios, and retains accessible table
headers. Its GitHub Actions run
https://github.com/EddyTodd/eddytodd-web/actions/runs/42172825404
passed Astro checks, all three benchmark-provenance tests, full build
and smoke. The website has **not** benchmarked untested later code.

## Pending DCT/FIR oracle and non-finite symmetry contract

The new `signal_processing_dct_fir_oracles` target independently
calculates unnormalized DCT-I/II/III/IV with textbook `long double`
cosine definitions (length zero and one, odd/prime/composite lengths
through 127, DC/impulse/alternating/oscillatory patterns). Its FIR
reference computes direct causal convolution with `long double`
accumulation for real and complex binary32/binary64, checking direct,
FFT, overlap-save, stateful streaming and reset responses across
short/empty and multi-block shapes.

An audit found the FIR symmetric/antisymmetric coefficient predicates
could incorrectly recognize NaN and infinity as valid due to unordered
floating-point comparisons and unchecked odd central taps; NaN/Inf
tolerances could also bypass a simple nonnegative guard. The predicates
now require **finite, nonnegative** tolerances and return false for
non-finite coefficients; the accelerated FIR paths consequently reject
those invalid coefficient arrays with `std::invalid_argument`.
The corresponding adversarial tests preserve the existing finite
coefficient numerical bounds. These changes have **not** passed new
Local CI yet: the next gate requires **12/12 Release** and **12/12
ASan/UBSan** executables and fresh reproducible timing artifacts.

## Second verified 12-suite campaign and correlation research extension

The latest pinned **same-source** report for
`589f871c8818f0db415e3b31f1d768659457a5e7` is
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/589f871c8818-afcb45f8-e5cb735bcbc66647
using coordinator `afcb45f8360ba5f11b8746d14f2563396f0197d3`. It passed
**12/12 Release** and **12/12 AddressSanitizer/UndefinedBehaviorSanitizer**
CTest executables, all with C assertions explicitly enabled, covering
the new long-double DCT/FIR regression suite and finite-input symmetry
contracts. The benchmark parser also passed. A full measured
**236-case** campaign was published (118 separately instantiated
binary32 + 118 binary64, each with nine raw timing samples for 2,124
total), i7-13700KF, GCC 13.3.0, `-O3 -march=native`, 102
CPU ISA capability flags, empirical median/p10/p90/MAD and matched
precision, shape and regime comparison rules. Artifact SHA-256:
`f72eba1712f8ef24c5da9e2bdd5b9a272f3ec7aecd259895b7cad30a4ade505c`.
The preceding report `589f871c8818-f4fee9dc-28e2506369435d61` used identical
source/benchmark build/compiler settings and passed the same gates,
but **all 236 timing medians differed** between machine runs. Timing
changes cannot establish an algorithm-level speedup or regression.

Website feature commit
`1a2ef1847f68e0c4bd6aca5517e8533b23e53c00` displays the
*previous* 12-suite artifact from `589f871c8818-f4fee9dc-28e2506369435d61`,
with genuine data, precision filtering, matched speedups and accessible
row-header semantics. It passed Astro checks (0 errors, 0 warnings),
three provenance tests, site build and smoke in
https://github.com/EddyTodd/eddytodd-web/actions/runs/42174428215.
The later same-source performance run is documented separately;
there is no aggregation, invented measurement or unsupported speedup
claim.

### Pending per-lag energy stability and convolution oracle

Audit found per-lag normalized correlation calculated overlap energy
from the difference of two large scalar running prefix sums. For a
single tiny trailing sample after a large sample, the difference can
round to zero even though the overlap has positive energy; the result
is incorrectly normalized to zero. A new `EnergyTree` uses pairwise
long-double summation *only* over the requested overlap range.
For FFT-based normalized correlation, the raw FFT cross-spectrum has
an absolute roundoff floor that can overwhelm a very weak overlap;
the implementation now adaptively recomputes affected lag numerators
by direct long-double products. The fallback is only for ill-conditioned
overlaps, and its threshold scales with the input energy and format
epsilon. It preserves the old non-FFT direct catalog API.

The new `signal_processing_convolution_oracles` independent CTest
uses long-double mathematical convolution/correlation references for
real and complex binary32/binary64 across empty, impulse, short, odd,
multi-block and sparse-tail inputs. It covers direct/FFT/circular
related paths, overlap-add, overlap-save, partitioned and streaming
flush/reset; invalid block sizes and weak normalized correlation
edges have explicit assertions and specimen diagnostics.
**These changes are unverified**, and the next laptop CI should
run **13/13 Release** and **13/13 ASan/UBSan** CTest executables,
then the full measured binary32/binary64 benchmark campaign.
No global mathematical guarantee, external performance baseline,
or all-input exceptional-value contract is claimed.
