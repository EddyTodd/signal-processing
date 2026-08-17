#pragma once

#include "signal_processing/fft.hpp"

#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace signal_processing::hilbert {

// Analytic signal / Hilbert transform via the DFT-domain multiplier.
// https://en.wikipedia.org/wiki/Analytic_signal
// https://en.wikipedia.org/wiki/Hilbert_transform

template <fft::Scalar T>
[[nodiscard]] inline std::vector<fft::Complex<T>> analytic_bluestein(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) return {};

    std::vector<fft::Complex<T>> spectrum(n);
    for (std::size_t i = 0; i < n; ++i) spectrum[i] = {input[i], T{}};
    spectrum = fft::bluestein<T>(spectrum, fft::Direction::forward);

    const std::size_t positive_end = (n + 1) / 2;
    for (std::size_t k = 1; k < positive_end; ++k) spectrum[k] *= T{2};
    const std::size_t negative_begin = n / 2 + 1;
    for (std::size_t k = negative_begin; k < n; ++k) spectrum[k] = {};
    // DC remains unchanged. For even N the Nyquist bin n/2 also remains unchanged.

    return fft::bluestein<T>(spectrum, fft::Direction::inverse);
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> transform_bluestein(std::span<const T> input) {
    const auto analytic = analytic_bluestein<T>(input);
    std::vector<T> output(analytic.size());
    for (std::size_t i = 0; i < analytic.size(); ++i) output[i] = analytic[i].imag();
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> envelope_bluestein(std::span<const T> input) {
    const auto analytic = analytic_bluestein<T>(input);
    std::vector<T> output(analytic.size());
    for (std::size_t i = 0; i < analytic.size(); ++i) output[i] = std::abs(analytic[i]);
    return output;
}

}  // namespace signal_processing::hilbert
