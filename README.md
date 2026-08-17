# signal-processing

A C++23 library of faithful implementations of important signal-processing algorithms. The repository is designed for both direct reuse and algorithm study: public APIs name the method being executed, and benchmarking infrastructure lives elsewhere.

## Scope

The current catalog contains:

- Fourier transforms: direct DFT, iterative/recursive radix-2, Stockham, radix-4, classical and modified split-radix, mixed-radix Cooley-Tukey, Good-Thomas/PFA, Rader, and Bluestein;
- reusable radix-2, real radix-2, mixed-radix, Good-Thomas, Rader, and Bluestein plans;
- radix 2/3/4/5/7 reusable codelets and explicit binary64 scalar/AVX2/AVX-512 radix-2 kernels;
- direct, circular, and FFT-based convolution;
- direct cross-correlation and autocorrelation;
- DCT-I, DCT-II, DCT-III, and DCT-IV;
- direct and streaming FIR filtering;
- Direct Form I IIR filtering and transposed Direct Form II biquads;
- Hann, Hamming, and Blackman windows.

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
const auto filtered = signal_processing::convolution::direct<double>(signal, kernel);
```

There is no algorithm registry or hidden dispatcher. Choose `fft::dft`, `fft::stockham`, `fft::mixed_radix`, `fft::rader`, `convolution::direct`, and so on explicitly. Reusable FFT plans are also algorithm-specific rather than selected by an opaque planner.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

The library is header-only and the CMake surface is intentionally small. Tests are built by default only when this repository is the top-level project. `tests/fft.cpp` checks every FFT mechanism against the direct DFT for deterministic binary32/binary64 inputs and validates plans, real transforms, inverse reconstruction, codelets, and available SIMD kernels.

## Layout

```text
include/signal_processing/         public algorithms
include/signal_processing/detail/  implementation support for complex algorithm families
tests/                             focused correctness and cross-method tests
docs/mathematics.md                definitions, conventions, and numerical notes
```

## Design rules

- C++23 with deliberate binary32 and binary64 support.
- Preserve algorithm identity: optimize implementations without silently changing the method.
- Prefer direct free functions and small stateful algorithm objects over registries and framework abstractions.
- Keep measurement, benchmark statistics, campaigns, and reporting outside this repository.
- Make normalization, supported lengths, state, scratch requirements, and ISA availability explicit.
- Architecture-specific implementations fail explicitly when unavailable; they never fall back to another kernel.
- Avoid host-dependent runtime `long double`; FFT chirp phase reduction uses overflow-safe modular integer arithmetic instead.
- Document mathematical definitions, normalization conventions, domains, and important numerical behavior.

## Remaining v1 milestones

1. convolution and correlation algorithms;
2. discrete transforms and window functions;
3. FIR filtering and multirate processing;
4. IIR filtering and filter design;
5. time-frequency and spectral analysis;
6. final numerical-contract and API consistency pass.

## Mathematics

See [`docs/mathematics.md`](docs/mathematics.md).

## License

MIT.
