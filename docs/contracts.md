# v1 numerical and API contracts

Version 1.0 freezes the public conventions below. Algorithm-specific mathematics remains in the other documents; this file defines the cross-library rules that callers can rely on.

## Scalar and sample formats

Public scalar numerical APIs use `float` or `double`. Host-dependent runtime `long double` is intentionally excluded.

Where the mathematics supports complex-valued samples, the corresponding sample types are `std::complex<float>` and `std::complex<double>`. Real-only algorithms, such as real linear-phase FIR design and IIR coefficient design, expose only the scalar formats.

Input **samples** follow ordinary IEEE floating-point propagation. A NaN or infinity in signal data is not silently replaced or rejected merely because it is nonfinite. In contrast, nonfinite values used to define algorithm geometry or coefficients—window parameters, filter coefficients, design ripple/attenuation, sample rates, normalized frequencies, and similar controls—are rejected when they do not define a meaningful algorithm.

## Explicit algorithm identity

Public method names identify the mechanism executed. There is no general dispatcher, implementation registry, hidden crossover threshold, or automatic benchmark-based selection.

Examples include `fft::stockham`, `fft::good_thomas`, `convolution::overlap_save`, `dct::type2_bluestein`, `fir::symmetric_direct`, `resampling::polyphase_rational`, and `spectral::welch_bluestein`.

Architecture-specific FFT kernels are also explicit. Requesting an unavailable AVX2 or AVX-512 implementation fails rather than falling back to scalar code.

## Errors

The library uses exceptions to make contract failures visible:

- `std::invalid_argument`: invalid domains, shapes, factors, coefficients, enum values, or incompatible buffers;
- `std::out_of_range`: a requested discrete index does not exist, such as a Goertzel bin outside the DFT;
- `std::length_error`: a required logical result, table, or workspace size is not representable in `std::size_t`;
- `std::runtime_error`: an iterative numerical method or structural numerical conversion fails to converge or satisfy its required invariant.

Algorithms do not intentionally continue with wrapped size arithmetic, unconverged exchange iterations, invalid SIMD substitutions, or malformed metadata.

## Empty inputs

Empty finite sequences are accepted where the mathematical operation has a natural empty result. Examples include DFT-family transforms, linear convolution with an empty operand, correlation with an empty operand, and STFT analysis of an empty signal.

An empty object is rejected when it would mean that the algorithm itself has not been defined. FIR execution therefore requires at least one coefficient, stateful convolution plans require a nonempty kernel, spectral density requires a nonempty positive-energy window, and a Goertzel **bin** requires an actual DFT bin. The arbitrary-frequency Goertzel sum of an empty sequence remains the mathematically valid zero sum.

## Size arithmetic

All size-sensitive algorithms must detect unrepresentable result/workspace sizes before wrapped arithmetic is used as an index or allocation size. This includes linear and blocked convolution, FFT convolution workspaces, Good-Thomas CRT maps and scratch storage, mixed-radix direct-leaf tables, multirate expansion, Remez grids, transformed IIR root counts, MDCT output sizing, and STFT frame offsets.

The implementation may still fail with the platform allocator for a representable but physically unavailable allocation; v1's contract is that the library itself does not intentionally derive that request through unsigned wraparound.

## Transform normalization

The complex DFT convention is

\[
X_k=\sum_{n=0}^{N-1}x_n e^{-2\pi i kn/N}.
\]

Forward transforms are unnormalized. Inverse FFT-family transforms apply exactly one factor of \(1/N\).

DCT, DST, Hartley, Walsh-Hadamard, and MDCT/IMDCT use the explicit unnormalized conventions in [`transforms.md`](transforms.md); their inverse scale factors are part of those definitions rather than a global orthonormal mode.

Good-Thomas follows the same call ordering as the other FFT methods:

```cpp
good_thomas(input, direction, factor_a, factor_b)
```

The factors may be omitted to use the algorithm's coprime factor split, but no different FFT family is selected.

## Frequency units

Frequency units are deliberately explicit even though different classical APIs conventionally use different normalizations:

- FIR design uses **cycles/sample**, with Nyquist equal to `0.5`;
- Goertzel arbitrary frequency also uses cycles/sample in `[0, 1)`;
- IIR digital design and IIR digital response use **normalized-to-Nyquist** frequency, with Nyquist equal to `1`;
- spectral-density estimators take an explicit sample rate and return density per sample-rate unit.

These units are not silently converted between APIs.

## FIR and multirate conventions

`fir::direct`, `fir::fft`, and `fir::overlap_save` return the causal response with one output per supplied input sample. `fir::full_direct` names the complete finite convolution explicitly.

Every FIR execution/resampling method requires a nonempty coefficient sequence. Symmetric and antisymmetric specialized implementations verify their coefficient structure and fail instead of falling back to the ordinary direct sum.

Rational/polyphase resampling does not silently reduce factors, invent a filter, scale interpolation gain, or choose a different implementation.

## IIR conventions

IIR coefficient arrays and explicit first/second-order sections require finite coefficients and nonzero `a0`. Constructors normalize by `a0` once.

Filter-design enum values are validated explicitly. High-order transfer-function coefficients remain available, while SOS conversion/execution is a separately named representation rather than a hidden substitution.

## Time-frequency and spectral conventions

STFT has no hidden centering or leading padding. `pad_end=true` emits only the minimal regular-hop final frame needed to cover the finite input. ISTFT performs weighted overlap-add and divides by the actual accumulated window-square weight.

Raw STFT power and spectral density are different APIs. PSD/CSD density scaling is

\[
\frac{1}{F_s\sum_n w[n]^2}.
\]

One-sided density is defined here only for real input; interior positive-frequency bins are doubled while DC and an even-length Nyquist bin are not.

Welch and CSD use complete frames only and do not silently detrend, subtract means, pad, or choose a window.

## Numerical equality

Faithful implementations of the same mathematical operation need not be bit-identical. FFT staging, blocking, partitioning, fused multiply-add, polyphase ordering, and other faithful reorganizations alter floating-point rounding order.

Correctness tests therefore compare against direct definitions or algebraic identities with binary32/binary64-specific tolerances. Exact bit equality is required only where an API or integer/index transformation specifically promises it.
