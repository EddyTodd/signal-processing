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

## Verified 13-suite normalized correlation and convolution campaign

The latest independently pinned `51c9e29045ce-91530d70-58022aa4ada6e37e` for source
`51c9e29045ce5891b77e0dc7dc1f9454a03095a0`, coordinator
`91530d7083be5ce2dd819f40295518494fbd1063`, finished
2026-11-06 06:45 UTC:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/51c9e29045ce-91530d70-58022aa4ada6e37e

All **13/13 Release** CTest executables passed with assertions active,
including the new `signal_processing_convolution_oracles` suite with
direct textbook `long double` convolution, real/complex binary32/binary64
circular/linear convolution, FFT implementations, overlap-add/save,
partitioned streaming and reset/flush, normalized correlation and weak-tail
regression. All **13/13 ASan/UBSan** CTest executables also passed;
the external benchmark parser self-tests passed.

The measured artifact has **236 cases** across FFT, convolution, DCT-II
and FIR (118 each for natively compiled binary32 and binary64), **nine
positive samples per case = 2,124 raw samples**, empirical median/p10/p90/MAD
and matched-regime speedups. CPU: 13th Gen Intel Core i7-13700KF,
compiler: GCC 13.3.0; flags `-O3 -DNDEBUG -std=c++23 -march=native
-fno-omit-frame-pointer`, 102 reported available CPU ISA flags (not a
disassembly/actual instruction-use audit). Published artifact SHA-256:
`df1caaa0e403a4890f96062560c0f7b29de6996241e42a781cef381e236c8d9b`.
Every case's nine-sample median, interpolated p10/p90, MAD, iteration
count, and matched baseline ratio was rechecked against the GitHub artifact.
No anomalous statistics or unmatched precision/regime speedup found.
The earlier `51c9e29045ce-e1d417be-ee89ec682ae2af39` has the same source, compiler flags, benchmark
code, machine and 13-suite PASS status but **all 236 median timings differ**.
Variation is not evidence of faster/slower algorithm changes.

The public website's feature commit
`95937f6a55c45453b2795b7e7e95537a99838fa5` pins the
earlier verified 13-suite artifact, not merged or invented timing data.
GitHub Actions
https://github.com/EddyTodd/eddytodd-web/actions/runs/42176097144
passed Astro checks, 3/3 benchmark provenance tests, full build,
and static/proxied route smoke tests. The newer same-source run is a
separate independent measurement, not an update of the displayed data.

## Pending finite IIR coefficient guard and fourth-order SOS oracle

IIR scalar normalization and first-order/biquad/SOS constructors
previously accepted `NaN`, `Inf`, or nonzero but sufficiently tiny
denominators whose quotient overflowed during normalization, silently
feeding invalid coefficients into stateful recursive filters. New
construction checks reject non-finite original or normalized
coefficients with `std::invalid_argument`, including nonfinite SOS
gain, while retaining legitimate tiny finite denominators when all
normalized ratios remain finite. This does **not** constrain sample
inputs or claim every finite filter is BIBO-stable.

New `signal_processing_iir_domains` checks invalid domains across
four direct/transposed realizations, first-order and three biquads,
and SOS constructors. A separately constructed fourth-order
two-section cascade is checked against long-double textbook recurrences
for impulse, step and oscillatory inputs, radius 0.76 and 0.985,
including reset and manually chained realizations. It extends,
rather than weakens, the previous IIR tests. **Not yet tested in
Local CI:** next report must show **14/14 Release and 14/14 ASan/UBSan**
with enabled assertions plus a new 236-case campaign. Neither
independent multi-machine performance evidence nor global
floating-point accuracy is established.

## Verified finite IIR and fourth-order SOS baseline (2026-11-06)

The pinned source commit `a8db0cbbb78fcef4cdba30457abc54696a44b9f7` and workflow
`1873d453ccf98ce9e38582ea70b25c303c25b9da` passed all stages in run `a8db0cbbb78f-1873d453-6ae773968c56a254`:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/a8db0cbbb78f-1873d453-6ae773968c56a254

**14/14 Release CTests** and **14/14 ASan/UBSan CTests**
passed with explicit `-UNDEBUG` assertions, including the
`signal_processing_iir_domains` invalid-coefficient and normalized
overflow contract for general, first-order, three biquad realizations and
SOS, and independent two-section fourth-order `long double` recurrence
against the rounded float32/float64 public coefficients. The benchmark
parser self-tests passed. The complete empirical benchmark contains
**236 cases** (118 float32 + 118 float64) and **2,124 timed samples**,
nine per case, on the Intel i7-13700KF with GCC 13.3.0 and
`-O3 -DNDEBUG -std=c++23 -march=native -fno-omit-frame-pointer`.
The reported 102 CPU ISA flags indicate *availability*, not actual
per-kernel instruction use. Every case's timing median, interpolated
p10/p90, median absolute deviation, trial count and same-shape/
same-precision/same-regime speedup was checked against the real samples.
Published artifact SHA-256:
`68f30f984b683c8cf22031f18cf592fcbdc789568b21d361ac0865e7ccad43bd`.

The eddytodd-web feature branch commit
`bc344e8475a96503d7dc4678805280190ff1f322`
pins that exact artifact, with accessible charts, filterable precision,
tabular data and measurement limits. GitHub Actions
https://github.com/EddyTodd/eddytodd-web/actions/runs/42222522016
passed Astro static analysis (0 errors/warnings), three data-provenance
tests, a full build and 23 static + one proxied route smoke validation.

### Pending ZPK and SOS real-coefficient input integrity

The public ZPK representation also permits callers to supply roots and
gain directly. An audit found that `second_order_sections` paired roots
within an approximate distance but then **discarded the imaginary part**
of the resulting quadratic's sum and product without checking its size.
Two non-conjugate roots close enough to match could therefore become
an apparently real biquad encoding a different transfer function.
Separately, non-finite roots/gain in transfer/SOS conversion or
non-finite analog evaluation frequency could enter calculations and
return non-finite values silently.

New input validation rejects NaN/Inf ZPK coordinates and gain and
non-finite analog evaluation frequencies, and rejects non-finite
expanded transfer coefficients and SOS section coefficients. SOS now
checks imaginary polynomial residuals against the existing
`4096 * epsilon` realness criterion before projecting to real
coefficients; this does not loosen the established polynomial
threshold. The extended `signal_processing_iir_domains` test checks
finite exact conjugate pairs, malformed approximate conjugate pairs,
NaN/Inf roots and gain, invalid analog frequencies and analog
stability with NaN imaginary poles for both float32 and float64.
This change is **not yet validated by laptop CI**: repeat all
**14/14 Release and 14/14 ASan/UBSan** executables and the full
precision-separated benchmark campaign. No full-space numerical
proof or external competitive performance comparison is claimed.

## Verified ZPK/SOS integrity, independent extreme-amplitude correlation regression pending

Pinned report `e53a2cefa111-1873d453-c0686d8e9a32d4cc` (tested source
`e53a2cefa111b1848b76a67050541ff2c1c56e14`, workflow `1873d453ccf98ce9e38582ea70b25c303c25b9da`) is
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/e53a2cefa111-1873d453-c0686d8e9a32d4cc.
All **14/14 Release** and **14/14 ASan/UBSan** CTests passed
with assertions active, including real-coefficient ZPK/SOS conjugate
integrity checks, non-finite ZPK roots/gains and analog frequencies,
fourth-order IIR recurrence, and all earlier independent transform,
spectral, filtering and convolution oracles. Benchmark-parser self-tests
passed.

Pinned benchmark SHA-256
`821f05b6d6eae4ee273bb4f607f57ec30062e58d1ba2ba1e73b86c1a47a8f9b1` contains **236 measured cases**
(118 separately compiled binary32 and 118 binary64), **2,124 actual
timed samples** (nine per method/shape/precision), correct medians,
interpolated p10/p90, MAD, and strictly same-precision, shape and regime
relative timing ratios. CPU Intel Core i7-13700KF, compiler GCC 13.3.0,
flags `-O3 -DNDEBUG -std=c++23 -march=native -fno-omit-frame-pointer`,
102 CPU ISA *available* flags (not proven executed ISA). CPU frequency,
cache, thermal state, external competitor, compiler portability and
multi-machine uncertainty are not controlled.

Website feature commit
`9eede12911dfaa2592991144b642992c3d935c5c`
ingests **exactly that measured artifact** with precision-selectable
accessible charts/tables, and GitHub workflow
https://github.com/EddyTodd/eddytodd-web/actions/runs/42223987117
passed Astro checks with zero errors/warnings, all three data-provenance
tests, complete build and 23 static plus one proxied route smoke tests.

### Pending normalized correlation dynamic-range repair

A separate audit found that normalizing only after forming the raw
cross-correlation in binary32/64 can give **incorrect zero or infinity**
when finite input amplitudes are too small/large for their product but
the normalized result is finite (e.g. self-correlation of a tiny
nonzero signal must be one at zero lag). The existing long-double
overlap-energy tree and FFT cancellation-sensitive fallback do not
automatically repair the **direct** path's intermediate underflow.
The new change recomputes an affected lag's numerator directly in long
double whenever the raw scalar result is non-finite or the
long-double normalization denominator lies outside the public scalar's
normal representable range. Independent textbook long-double oracles
now test extreme finite amplitudes for real and complex float32/float64,
including conjugation and far-lag single-tap overlaps. The test bounds
remain unchanged for earlier cases; the new adversarial cases have
explicit independent expected values. **This change is not yet tested
on the laptop**; the next gate must verify all **14/14 Release and
14/14 ASan/UBSan CTests**, then a fresh benchmark/parser campaign.
The range of `long double` differs by platform; this is not a
portable guarantee that all double extremes normalize accurately.

## Verified extreme-range normalized correlation and next spectral domain gate

Pinned report `d07d6db4f24c-1873d453-0b30aba5b2261519` for source `d07d6db4f24cd4c111bc466c39c17457dc74c6b4` and coordinator
`1873d453ccf98ce9e38582ea70b25c303c25b9da` passed all stages:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/d07d6db4f24c-1873d453-0b30aba5b2261519

**14/14 Release** and **14/14 ASan/UBSan** CTest executables passed with
`-UNDEBUG` assertions enabled, including direct and FFT normalized
correlation of very small/large real/complex binary32 and binary64 values,
plus independent long-double mathematical oracles and all earlier
IIR/ZPK, DCT, FIR, FFT, convolution and STFT/spectral suites. The
benchmark parser self-tests passed. Benchmark artifact SHA-256:
`f87d3953a81b0e42c472e704b5faaa6416cad0f59cc899ab9fb2b28a4fb8fce2`: **236 measured cases** (118 float32, 118 float64), **2,124
actual raw timing samples** (nine per case), independent checks of
median/p10/p90/MAD and matched speedup ratios. CPU i7-13700KF and
GCC 13.3.0, flags `-O3 -DNDEBUG -std=c++23 -march=native
-fno-omit-frame-pointer`; available ISA flags recorded, instruction use
unproven. Timings are noisy empirical observations, not confidence
intervals or cross-machine performance claims.

Website feature commit `5b46e70a47664e4cbb0e452ef7f60ab2c1f5f266`
pins only this exact measured artifact and its test provenance. GitHub
Actions run https://github.com/EddyTodd/eddytodd-web/actions/runs/42225998150
passed Astro checks (0 errors and warnings), all three benchmark
provenance tests, full build and 23-static/1-proxied route smoke.

### Pending checked spectral parameter/normalization contracts

Audit found that the previous PSD implementations accepted `+Inf`
sample rate because it satisfied the positive-only comparison, and
could silently output zeros when `sample_rate * window_energy`
overflowed. Non-finite window entries could also produce non-finite
spectral density. An invalid `Sides` enumerator was treated like
two-sided, while a mismatched window on an empty periodogram or a
malformed window on a too-short Welch/CSD signal could bypass validation.

The new `spectral.hpp` guards require finite, positive sample rates,
supported side values, finite window coefficients and positive finite
window energy. The normalization product and reciprocal scale must
both be finite and strictly positive in the public binary32/64
format; otherwise they explicitly throw `std::invalid_argument`
rather than return silently corrupted density. Empty matching
periodogram inputs continue to return empty; a shape mismatch
throws, and short Welch/CSD still return empty *after* validating
window energy. New `signal_processing_spectral_domains` tests four
real/complex sample types, invalid domain and normalization extremes
and early-return boundaries. **This new code and 15th CTest target
are not yet validated on the laptop:** the next report must show
**15/15 Release** and **15/15 ASan/UBSan** with assertions active,
the benchmark/parser campaign and source-SHA-matched artifact.
Finite input/output guardrails do not prove all intermediate FFT
products remain representable; no full numerical stability or
third-party competitiveness is claimed.

## Verified 15-suite spectral domain campaign (2026-11-06)

Pinned source `0ee79b23ee013e0131b32d0decfea5639d82a0de`, workflow `a5ca79128a076396583962b39b22d4db6b5321ce`, run `0ee79b23ee01-a5ca7912-9aeb772c1efca341`:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/0ee79b23ee01-a5ca7912-9aeb772c1efca341

All **15/15 Release** and **15/15 ASan/UBSan** CTest executables
passed with `-UNDEBUG` assertions enabled, including independent
real/complex float32/float64 spectral-domain checks for NaN/Inf/invalid
rate and window inputs, empty/short input contracts, unrepresentable
spectral density normalization and unknown sides enum. Earlier
independent long-double Fourier, IIR/ZPK, STFT, convolution, DCT/FIR
regressions also passed. Benchmark parser self-tests passed.

The pinned artifact (SHA-256 `195274472d4b9332964b51e899038022ba6d5d780d6a3b3f8e95593a7f99a1f1`) includes **236
actual method/shape cases**, 118 binary32 and 118 binary64, with **9
raw timed samples per case = 2,124 total**, verified medians, p10/p90,
MAD and matched precision/shape/regime baseline ratios. GCC 13.3.0
and Intel i7-13700KF, `-O3 -DNDEBUG -std=c++23 -march=native
-fno-omit-frame-pointer`, 102 *available* ISA flags (not a per-kernel
disassembly audit). Spread is empirical single-machine variability,
not a confidence interval. No external FFTW/SciPy, other-compiler,
other-ISA or thermal-controlled performance results are implied.

Website feature commit `dbae6a84c6dbc786f4398ea6c928bd84bd724d29` displays **exactly this**
GitHub-published source-matched artifact and accessible
precision-selectable comparisons. GitHub Actions
https://github.com/EddyTodd/eddytodd-web/actions/runs/42227671594
PASSED Astro (0 errors/warnings), all 3/3 data-provenance tests,
complete build and 23 static + 1 proxied smoke routes.

### Pending STFT inverse weight and overflow contracts

A direct audit identified acceptance of nonfinite STFT analysis-window
weights even for empty inputs, plus inverse analysis on user-built
frames with NaN/Inf bins or on windows whose squared weights underflow
or overflow. Even if each weight-square remains finite, many overlapping
frames can overflow their sum and silently poison inverse output.

The new `stft.hpp` guards validate window coefficients and inverse
frame-bin finiteness, reject squared nonzero weights that collapse to
zero or infinity, and detect nonfinite accumulated inverse overlap
numerators or denominators. Intentional zero windows preserve zero
reconstruction. A new `signal_processing_stft_domains` 16th CTest
covers these exceptional domains and valid finite negative-window
reconstruction across real/complex binary32 and binary64.
**Not yet verified on laptop**: the next run must pass **16/16 Release**
and **16/16 ASan/UBSan**, then rerun and publish 236 full dual-precision
nine-sample benchmark cases plus parser self-tests. Remaining gaps
include pathological intermediate FFT overflow, near-unit high-order
IIR conditioning, complete API-wide NaN/Inf behavior, external
benchmark baselines and cross-machine reproducibility.

## Two verified 16-suite STFT domain campaigns and pending analytic oracles

Source `584a2878254ea4b2b287be19153ab5f7bb15af35` passed **16/16 Release** and **16/16 ASan/UBSan**
CTest executables with `-UNDEBUG` assertions active in two distinct
pinned Local CI reports:

- https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/584a2878254e-594b7546-54d69b3205978627
- https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/584a2878254e-b9344b34-edce10bd802439ad

The newer workflow configuration `b9344b3438397b8f7afc500d6b61cafc73f127c4` used the same benchmark
source code, GCC 13.3.0 and Intel i7-13700KF, and again passed all
four stages including benchmark-parser self-tests. Measured campaign
artifact SHA-256 `ad2f3224085072652d80f1c9aec60d8650c21c02582baa416fbd48d6acc82a63` contains **236** native binary32/binary64
FFT/convolution/DCT/FIR cases (118 shapes per precision), **nine
timings per case** = **2,124 original samples**. Median, interpolated
p10/p90, MAD and matched-precision/shape/regime speedups were
independently checked against raw published samples. All 236 timing
medians differ from the earlier same-source report even though the
benchmark kernels and compiler flags were identical, so they are not
evidence of code-level performance change. Available CPU ISA flags
do not establish actual instructions executed. Results are empirical
single-machine observations, not a cross-platform confidence interval.

Website feature commit `060c162b17f6df82722606ead7a83559a193b466` deliberately retains the
earlier verified, measured 16-suite source-matched artifact
`5ae3e69f7fafd7aa7239fe380f21a9973a0ae761ba1255b3f5bd155ef8e85763`,
without inventing an aggregate across noisy runs.
https://github.com/EddyTodd/eddytodd-web/actions/runs/42230381911
succeeded: Astro 0 errors/0 warnings, 3/3 provenance tests, full
build, and 23 static + 1 proxied route smoke tests.

### Next independent analytic-signal numerical gate

Existing Goertzel and Hilbert tests covered relatively few lengths,
off-bin frequencies or odd-length analytic signals. An additional
`signal_processing_analytic_oracles` executable was committed and
registered as the **17th** CTest. It constructs a mathematically
independent `long double` DFT sum for Goertzel coefficients and an
independent direct Fourier-domain Hilbert analytic-signal transform,
covering complex and real float32/float64 signals, impulse,
oscillatory/chirp-like inputs, empty/odd/even/prime/composite lengths,
DFT bins, nonbin frequencies near Nyquist and 1, and signal envelopes.

The new test is **not** part of either passed 16-suite run; next laptop
Local CI must compile and pass **17/17 Release and 17/17 ASan/UBSan**
with assertions enabled and publish the complete 236-case dual-precision
benchmark/parser campaign. External FFTW/SciPy baselines,
cross-compiler/ISA/platform replication, extreme IIR conditioning,
global numerical proof and controlled machine uncertainty remain open.

## Repeated validated seventeen-suite analytic-oracle campaign

Source `68d5e401e63fbdce9392e5f316a1a13144d06668` has two independently pinned PASS reports:
- https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/68d5e401e63f-f5c623d8-f9d9c92f52a5f0c4
- https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/68d5e401e63f-d2e69a43-72b81618b891aad9

The newer configuration `d2e69a4341c2034acf1d10f5afe03aa0c4d71ed2` passed **17/17 Release** and
**17/17 ASan/UBSan** CTest executables with assertions explicitly active
(`-UNDEBUG`), including independent direct `long double` Fourier
Goertzel complex phase and Hilbert analytic/envelope odd/even
regressions, all earlier spectral/FFT/IIR/FIR/convolution/DCT oracles,
and benchmark parser self-tests.

The published second artifact SHA-256 `fc70a02845dfc9b6427ebfda4014f7e555d7d103e5576a1b548eec6fe35989cc` contains **236**
native-precision FFT/convolution/DCT/FIR method-shapes, 118 binary32
and 118 binary64, **2,124 timed samples** (nine each), with all reported
median, interpolated p10/p90, MAD and matched precision/shape/regime
ratios independently checked. On the same i7-13700KF laptop,
GCC 13.3.0, `-O3 -DNDEBUG -std=c++23 -march=native
-fno-omit-frame-pointer`, with 102 available CPU ISA flags,
**all 236 medians differ** from the earlier same-source, same-kernel
run. This is run-to-run timing variability; it cannot establish a
speedup or regression. Available ISA metadata does not verify which
instructions executed.

Website branch `feature/signal-processing-lab-20261105` commit
`525486f2954cea65f881ebe5e4c68e00c3f59762` pins the
*earlier* verified 17-suite source-matched artifact SHA-256
`af716476c494388f9a6a7add2de4028b35d773bb7de987fc6680811aedad091e`.
It is not merged or averaged with the newer measurement.
GitHub workflow
https://github.com/EddyTodd/eddytodd-web/actions/runs/42236331082
passed Astro 0 errors/0 warnings, all three benchmark provenance
tests, full build and 23 static plus 1 proxied smoke routes.

### Pending Goertzel/Hilbert exceptional-value contract

Even after independent phase and analytic-signal oracles passed, the
public Goertzel recurrence accepted NaN/Inf real or imaginary input
samples and could return nonfinite states or coefficients from huge
but finite input. Hilbert analytic, transform and envelope calls also
accepted nonfinite real samples. New source checks reject nonfinite
Goertzel inputs with `std::invalid_argument`, recurrence/final
coefficient overflow with `std::overflow_error`, and nonfinite
Hilbert real input and resulting analytic complex samples. Existing
`signal_processing_analytic_oracles` additionally checks both public
Goertzel entry points, complex imaginary components, known finite
overflow, and all three Hilbert entry points across binary32/binary64.
The previous independent finite-oracle tolerances were not relaxed.
**This newer source has not passed laptop CI.** The next report must
verify **17/17 Release and 17/17 ASan/UBSan**, assertions active,
and the full measured dual-precision benchmark/parser suite.

There is no proof of bounded error for every finite coefficient,
no full API-wide exceptional-value certification or high-order IIR
conditioning guarantee. External FFTW/SciPy comparisons, alternate
compilers/ISA/platforms, controlled thermal/frequency/caches and
independent multi-machine performance remain outstanding.

## Verified finite analytic-contract campaign; pending spectral-data validation

The exact pinned CI run `92ca02aed38a-0de4ef4b-f3d9589fbd740635` for source
`92ca02aed38ac6f1df73ba92ef423aac455caf6b` and coordinator `0de4ef4b4209fbbe41c3d27cbf867ad7e78b2818` passed all
**17/17 Release** and **17/17 ASan/UBSan** CTest executables
with assertions explicitly active (`-UNDEBUG`), including independent
Goertzel/Hilbert long-double phase and analytic checks plus nonfinite
real/complex sample rejection and overflowing recurrence/analytic-result
checks. All four stages, including benchmark parser self-tests, passed:
https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/92ca02aed38a-0de4ef4b-f3d9589fbd740635

The GitHub benchmark artifact SHA-256 `3e6fd4fe9efd373223a486f26fabf362761863ac9e5637e6818c72f9afaeddd8` contains
**236 real FFT/convolution/DCT/FIR method-shape benchmarks** (118
binary32 and 118 binary64), **2,124 actual timings** (nine per case).
Every reported median, p10/p90, MAD, and matched
precision/shape/allocation-regime baseline ratio was independently
checked against published raw samples. CPU: Intel i7-13700KF, GCC 13.3.0,
`-O3 -DNDEBUG -std=c++23 -march=native -fno-omit-frame-pointer`.
102 ISA capabilities were reported; executed ISA was not audited.
Observed spread is empirical and does not establish an algorithmic
speedup or cross-platform confidence interval.

The eddytodd-web feature commit
`b0a5e512c04398510d0122a37a7b12acba932b02` displays exactly
that source-matched artifact. GitHub Actions
https://github.com/EddyTodd/eddytodd-web/actions/runs/42238308945
passed Astro (zero errors and warnings), all three benchmark-data
provenance tests, full build and 23 static + one proxied route smoke.

### Pending finite-spectral-sample and density-output integrity

Audit of the public periodogram, Welch, CSD and spectrogram estimators
found that they validated rates/windows but not real or complex sample
values. NaN/Inf inputs could become NaN/Inf FFT coefficients and
densities, while large **finite** inputs could overflow complex
multiplications, FFT bins, powers or one-sided PSD doubling. Malformed
samples in short Welch/CSD inputs could even be ignored on the early
empty-result path.

New guards validate every real/imaginary input sample, including short
or otherwise zero-frame input, check finite windowed FFT intermediates,
and reject nonfinite PSD/Welch/CSD/spectrogram outputs. The independently
expanded `signal_processing_spectral_domains` CTest checks those
contracts for native binary32/binary64 real/complex data, including
finite values too large for normalized Fourier or density arithmetic.
The new tests do not loosen prior numeric tolerances. **Not yet run
in Local CI**: next report must confirm **17/17 Release and 17/17
ASan/UBSan** with assertions enabled plus the complete 236-case
nine-sample-per-case native-precision benchmark and parser self-tests.

Remaining research limits: exceptional-value coverage throughout all
other DSP APIs, finite-but-inaccurate output, high-order near-unit
pole conditioning, external FFTW/SciPy baseline and alternative
compiler/ISA/platform replication, controlled CPU conditions, memory
throughput and independent machines. Do not mark COMPLETE.

## Two verified 17-suite spectral-input campaigns; pending FIR symmetry repair

Tested source `aaf3a39e506d777103a0da131829048b1cee2c91`, same pinned compiler/benchmarks, two
independent PASS reports:
- https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/aaf3a39e506d-0de4ef4b-87cbf3ccee1ac42f
- https://github.com/EddyTodd/local-ci-reports/tree/main/reports/autonomous/EddyTodd/signal-processing/ai-localci-dsp-research/aaf3a39e506d-7615bdc1-b4c6f3741e1efba2

The newer coordinator SHA `7615bdc1052efec41939a42294eb90614c7d7f6f` again passed all
**17/17 Release** and **17/17 ASan/UBSan** CTest executables with
`-UNDEBUG` assertions active, including spectral real/complex NaN/Inf
samples, finite-but-overflowing FFT/PSD intermediates, short Welch/CSD
and spectrogram input-domain guards, and all previous Goertzel/Hilbert,
STFT, FIR/DCT, IIR, FFT and convolution numerical oracles. Benchmark
parser self-tests passed. GitHub artifact SHA-256
`7e9ce4cf86fc0213cf5428e862bab54d1e4ddd309c95ca99cf44186af0d42f60`: **236 measured native-precision benchmarks** (118
binary32 + 118 binary64) and **2,124 original nine-trial samples**;
every median, p10, p90, MAD and precision/shape/regime-matched
baseline speedup was checked against the real samples. Intel i7-13700KF,
GCC 13.3.0, `-O3 -DNDEBUG -std=c++23 -march=native
-fno-omit-frame-pointer`, 102 reported available CPU ISA flags
(which do not establish instructions actually used). All 236 medians
differ from the prior same-source measurement; this is uncontrolled
run-to-run timing variability, not evidence of code-level acceleration.

The website feature commit `a0e0c820e924743bb7b7c5be386aa100968c6869` pins the *earlier*
independently verified `aaf3a39e506d-0de4ef4b-87cbf3ccee1ac42f` source-matched measured artifact
`89454c45da2af5f52e9f0138878d61794d20c78081abb1fa90a929d993ef0b5e`,
with accessible charts and tables, and does not blend the runs.
https://github.com/EddyTodd/eddytodd-web/actions/runs/42240292689
PASSED Astro checks (zero errors and warnings), three benchmark-data
provenance tests, complete build, and 23 static plus one proxied route
smoke.

### Pending overflow-safe linear-phase FIR coefficient classification

In `fir::symmetric_coefficients` and
`fir::antisymmetric_coefficients`, two finite but near-maximum taps
could overflow `std::abs(a-b)` or `std::abs(a+b)` **and**
the original `tolerance*(1+max(abs(a),abs(b)))` bound. Comparing
`Inf > Inf` incorrectly produced `false`, silently accepting an
invalid linear-phase tap vector and dispatching a mathematically
different symmetric/antisymmetric optimized filter.

The new scaled residual computation divides both finite taps by
`max(1,abs(a),abs(b))` before addition/subtraction, then compares
their normalized pair mismatch against the original tolerance
without constructing a potentially overflowing intermediate.
The independent regression checks max finite opposite/equal/unbalanced
pairs, negative taps, zero and very large tolerances, tiny subnormal
values, and rejection by both optimized FIR entry points for real
binary32/binary64. Earlier finite numerical-error bounds were not
changed. **This new source has not yet run laptop CI:** the next
report must confirm **17/17 Release** and **17/17 ASan/UBSan**, the
parser self-tests and the full 236-case native binary32/binary64
nine-timing campaign before updating website attribution.

Research limits: remaining API-wide exceptional-value cases,
near-unit high-order IIR conditioning, global numerical proof,
external FFTW/SciPy and cross-compiler/ISA/OS benchmarks, CPU
thermal/cache/frequency controls and multi-machine replication.
Keep task ACTIVE.
