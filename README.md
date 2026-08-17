# signal-processing

A C++23 library of faithful implementations of important signal-processing algorithms. The repository is designed for both direct reuse and algorithm study: public APIs name the method being executed, and benchmarking infrastructure lives elsewhere.

## Scope

The initial catalog contains:

- direct DFT and iterative radix-2 FFT;
- direct, circular, and FFT-based convolution;
- direct cross-correlation and autocorrelation;
- DCT-I, DCT-II, DCT-III, and DCT-IV;
- direct and streaming FIR filtering;
- Direct Form I IIR filtering and transposed Direct Form II biquads;
- Hann, Hamming, and Blackman windows.

The FFT catalog will be expanded from the algorithms already researched in [`EddyTodd/fft`](https://github.com/EddyTodd/fft), but FFT remains part of this repository rather than an external runtime dependency.

## Use

```cpp
#include <signal_processing/signal_processing.hpp>

#include <complex>
#include <vector>

std::vector<std::complex<double>> values{{1.0, 0.0}, {2.0, 0.0},
                                         {3.0, 0.0}, {4.0, 0.0}};
const auto spectrum = signal_processing::fft::radix2<double>(values);

std::vector<double> signal{1.0, 2.0, 3.0};
std::vector<double> kernel{0.5, 0.5};
const auto filtered = signal_processing::convolution::direct<double>(signal, kernel);
```

There is no algorithm registry or hidden dispatcher. Choose `fft::dft`, `fft::radix2`, `convolution::direct`, `convolution::fft`, and so on explicitly.

## Build and test

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

The library is header-only and the CMake surface is intentionally small. Tests are built by default only when this repository is the top-level project.

## Layout

```text
include/signal_processing/  public algorithms
tests/                      focused correctness tests
docs/mathematics.md         definitions, conventions, and numerical notes
```

## Design rules

- C++23 with explicit `float` and `double` support through constrained templates.
- Preserve algorithm identity: optimize implementations without silently changing the method.
- Prefer direct free functions and small stateful filter objects over registries and framework abstractions.
- Keep measurement, benchmark statistics, campaigns, and reporting outside this repository.
- Keep tests focused on correctness and cross-method agreement rather than exhaustive benchmark infrastructure.
- Document mathematical definitions, normalization conventions, domains, and important numerical behavior.

## Planned extensions

The next algorithm-focused milestones include the remaining representative FFT mechanisms, overlap-add and overlap-save convolution, FFT-based DCTs, additional IIR forms and second-order-section cascades, FIR optimizations, DSTs, resampling/polyphase methods, STFT, and spectral analysis.

## Mathematics

See [`docs/mathematics.md`](docs/mathematics.md).

## License

MIT.
