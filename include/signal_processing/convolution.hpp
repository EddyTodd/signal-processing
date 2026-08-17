#pragma once

#include "signal_processing/fft.hpp"

#include <algorithm>
#include <complex>
#include <concepts>
#include <cstddef>
#include <span>
#include <vector>

namespace signal_processing::convolution {

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> direct(std::span<const T> signal,
                                            std::span<const T> kernel) {
    if (signal.empty() || kernel.empty()) {
        return {};
    }
    std::vector<T> output(signal.size() + kernel.size() - 1, T{});
    for (std::size_t i = 0; i < signal.size(); ++i) {
        for (std::size_t j = 0; j < kernel.size(); ++j) {
            output[i + j] += signal[i] * kernel[j];
        }
    }
    return output;
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> circular(std::span<const T> lhs,
                                              std::span<const T> rhs) {
    if (lhs.empty() || lhs.size() != rhs.size()) {
        return {};
    }
    const std::size_t n = lhs.size();
    std::vector<T> output(n, T{});
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            output[k] += lhs[j] * rhs[(k + n - j) % n];
        }
    }
    return output;
}

[[nodiscard]] constexpr std::size_t next_power_of_two(std::size_t n) noexcept {
    if (n <= 1) {
        return 1;
    }
    --n;
    for (std::size_t shift = 1; shift < sizeof(std::size_t) * 8; shift <<= 1) {
        n |= n >> shift;
    }
    return n + 1;
}

// Convolution theorem.
// https://en.wikipedia.org/wiki/Convolution_theorem
template <std::floating_point T>
[[nodiscard]] inline std::vector<T> fft(std::span<const T> signal,
                                         std::span<const T> kernel) {
    if (signal.empty() || kernel.empty()) {
        return {};
    }
    const std::size_t output_size = signal.size() + kernel.size() - 1;
    const std::size_t n = next_power_of_two(output_size);
    using Complex = signal_processing::fft::Complex<T>;
    std::vector<Complex> a(n);
    std::vector<Complex> b(n);
    for (std::size_t i = 0; i < signal.size(); ++i) {
        a[i] = Complex{signal[i], T{0}};
    }
    for (std::size_t i = 0; i < kernel.size(); ++i) {
        b[i] = Complex{kernel[i], T{0}};
    }

    signal_processing::fft::radix2_inplace<T>(a);
    signal_processing::fft::radix2_inplace<T>(b);
    for (std::size_t i = 0; i < n; ++i) {
        a[i] *= b[i];
    }
    signal_processing::fft::radix2_inplace<T>(a, signal_processing::fft::Direction::inverse);

    std::vector<T> output(output_size);
    std::transform_n(a.begin(), output_size, output.begin(),
                     [](const Complex& value) { return value.real(); });
    return output;
}

}  // namespace signal_processing::convolution
