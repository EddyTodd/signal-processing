#pragma once

#include "signal_processing/fft.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <vector>

namespace signal_processing::hartley {

// Discrete Hartley transform.
// https://en.wikipedia.org/wiki/Discrete_Hartley_transform

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> direct(std::span<const T> input) {
    const std::size_t n = input.size();
    std::vector<T> output(n, T{});
    if (n == 0) return output;

    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k = 0; k < n; ++k) {
        T sum{};
        for (std::size_t j = 0; j < n; ++j) {
            const T angle = tau * static_cast<T>(j) * static_cast<T>(k) /
                            static_cast<T>(n);
            sum += input[j] * (std::cos(angle) + std::sin(angle));
        }
        output[k] = sum;
    }
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> bluestein(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) return {};

    std::vector<fft::Complex<T>> complex_input(n);
    for (std::size_t j = 0; j < n; ++j) complex_input[j] = {input[j], T{0}};
    const auto spectrum = fft::bluestein<T>(complex_input);

    std::vector<T> output(n);
    for (std::size_t k = 0; k < n; ++k) {
        output[k] = spectrum[k].real() - spectrum[k].imag();
    }
    return output;
}

}  // namespace signal_processing::hartley
