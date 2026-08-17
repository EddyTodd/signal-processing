# signal-processing

A C++23 library of faithful implementations of important signal-processing algorithms. The repository is designed for both direct reuse and algorithm study: public APIs name the method being executed, and benchmarking infrastructure lives elsewhere.

## Scope

The current catalog contains:

- Fourier transforms: direct DFT, iterative/recursive radix-2, Stockham, radix-4, classical and modified split-radix, mixed-radix Cooley-Tukey, Good-Thomas/PFA, Rader, and Bluestein;
- reusable radix-2, real radix-2, mixed-radix, Good-Thomas, Rader, and Bluestein plans;
- radix 2/3/4/5/7 reusable codelets and explicit binary64 scalar/AVX2/AVX-512 radix-2 kernels;
- real and complex direct, circular, radix-2 FFT, overlap-add, overlap-save, streaming direct, and uniform partitioned convolution;
- real and complex direct/FFT cross-correlation, autocorrelation, and per-lag normalized correlation;
- DCT-I through DCT-IV as direct definitions and arbitrary-length Bluestein-derived forms;
- DST-I through DST-IV as direct definitions and arbitrary-length Bluestein-derived forms;
- direct/Bluestein discrete Hartley transforms and direct/fast Sylvester-ordered Walsh-Hadamard transforms;
- direct MDCT and IMDCT definitions;
- rectangular, Bartlett, Hann, Hamming, Blackman, Blackman-Harris, Nuttall, flat-top, Kaiser, Gaussian, Tukey, Lanczos, and Welch windows with explicit symmetric/periodic sampling where applicable;
- real/complex direct, FFT, overlap-save, streaming, symmetric, and antisymmetric FIR filtering;
- windowed-sinc low/high/band-pass/band-stop FIR design and Type-I Parks-McClellan/Remez equiripple design;
- downsampling, zero insertion, zero-order hold, linear and windowed-sinc interpolation, FIR decimation/interpolation, rational conversion, and explicit polyphase forms;
- Direct Form I IIR filtering and transposed Direct Form II biquads.

The FFT implementation is owned by this repository. `EddyTodd/fft` remains useful as the earlier focused research history, but `signal-processing` does not depend on it.

## Use

```cpp
#include <signal_processing/signal_processing.hpp>

#include <complex>
#include <vector>

std::vector<std::complex<double>> values{{1.0, 0.0}, {2.0, 0.0},
                                         {3.0, 0.0}, {4.0, 0.0}};
const auto spectrum = signal_processing::fft::radix2<double>(values);

signal_processing::fft::Radix2Plan<double> plan(1024);
std::vector<std::complex<double>> block(1024);
plan.forward_inplace(block);

std::vector<double> signal{1.0, 2.0, 3.0};
std::vector<double> kernel{0.5, 0.5};
const auto direct = signal_processing::convolution::direct<double>(signal, kernel);
const auto blocked = signal_processing::convolution::overlap_add<double>(signal, kernel, 2);

const auto cosine = signal_processing::dct::type2<double>(signal);
const auto cosine_fast = signal_processing::dct::type2_bluestein<double>(signal);
const auto analysis_window = signal_processing::windows::hann<double>(
    1024, signal_processing::windows::Sampling::periodic);

const auto lowpass = signal_processing::fir_design::remez_lowpass<double>(
    63, 0.18, 0.24, 1.0, 10.0);
const auto converted = signal_processing::resampling::polyphase_rational<double>(
    signal, 3, 2, lowpass);
```

There is no algorithm registry or hidden dispatcher. Choose `fft::dft`, `fft::stockham`, `fft::mixed_radix`, `convolution::direct`, `convolution::overlap_add`, `dct::type2`, `dct::type2_bluestein`, `fir::direct`, `fir::overlap_save`, `resampling::rational_fir`, `resampling::polyphase_rational`, and so on explicitly. Reusable state is algorithm-specific rather than selected by an opaque planner.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

The library is header-only and the CMake surface is intentionally small. Tests are built by default only when this repository is the top-level project. `tests/fft.cpp` checks the FFT mechanisms against the direct DFT. `tests/convolution_correlation.cpp` checks blocked, streaming, real/complex, and frequency-domain convolution/correlation forms against direct definitions. `tests/transforms.cpp` checks direct/fast transform agreement, inverse-scaling identities, MDCT basis behavior, and window invariants. `tests/fir_multirate.cpp` checks direct/FFT/overlap-save/streaming FIR equivalence, linear-phase reductions, FIR design invariants, and direct/polyphase multirate equivalence.

## Layout

```text
include/signal_processing/         public algorithms
include/signal_processing/detail/  shared implementation support
tests/                             focused correctness and cross-method tests
docs/mathematics.md                core definitions and numerical conventions
docs/transforms.md                 transform/window formulas and normalization
docs/fir-multirate.md              FIR design/execution and sample-rate conversion
```

## Design rules

- C++23 with deliberate binary32 and binary64 support, including real and complex signal samples where the mathematics supports both.
- Preserve algorithm identity: optimize implementations without silently changing the method.
- Prefer direct free functions and small stateful algorithm objects over registries and framework abstractions.
- Keep measurement, benchmark statistics, campaigns, and reporting outside this repository.
- Keep design algorithms separate from execution algorithms.
- Make normalization, supported lengths, transform conventions, block/partition sizes, resampling factors, state, scratch requirements, and ISA availability explicit.
- Architecture-specific implementations fail explicitly when unavailable; they never fall back to another kernel.
- Avoid host-dependent runtime `long double`; FFT chirp phase reduction uses overflow-safe modular integer arithmetic instead.
- Document mathematical definitions, normalization conventions, domains, and important numerical behavior.

## Remaining v1 milestones

1. IIR filtering and filter design;
2. time-frequency and spectral analysis;
3. final numerical-contract and API consistency pass.

## Mathematics

See [`docs/mathematics.md`](docs/mathematics.md) for core algorithms, [`docs/transforms.md`](docs/transforms.md) for transform/window conventions, and [`docs/fir-multirate.md`](docs/fir-multirate.md) for FIR and multirate processing.

## License

MIT.
