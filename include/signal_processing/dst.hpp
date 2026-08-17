#pragma once

#include "signal_processing/fft.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_processing::dst {

// Discrete sine transform.
// https://en.wikipedia.org/wiki/Discrete_sine_transform

namespace detail {

[[nodiscard]] inline std::size_t checked_twice(std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / 2) {
        throw std::length_error("DST workspace size overflow");
    }
    return 2 * n;
}

}  // namespace detail

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> type1(std::span<const T> input) {
    const std::size_t n = input.size();
    std::vector<T> output(n, T{});
    for (std::size_t k = 0; k < n; ++k) {
        T sum{};
        for (std::size_t j = 0; j < n; ++j) {
            const T angle = std::numbers::pi_v<T> * (static_cast<T>(j) + T{1}) *
                            (static_cast<T>(k) + T{1}) / (static_cast<T>(n) + T{1});
            sum += input[j] * std::sin(angle);
        }
        output[k] = T{2} * sum;
    }
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> type2(std::span<const T> input) {
    const std::size_t n = input.size();
    std::vector<T> output(n, T{});
    if (n == 0) return output;

    for (std::size_t k = 0; k < n; ++k) {
        T sum{};
        for (std::size_t j = 0; j < n; ++j) {
            const T angle = std::numbers::pi_v<T> * (static_cast<T>(j) + T{0.5}) *
                            (static_cast<T>(k) + T{1}) / static_cast<T>(n);
            sum += input[j] * std::sin(angle);
        }
        output[k] = T{2} * sum;
    }
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> type3(std::span<const T> input) {
    const std::size_t n = input.size();
    std::vector<T> output(n, T{});
    if (n == 0) return output;

    for (std::size_t k = 0; k < n; ++k) {
        T sum = (k & 1U) == 0U ? input.back() : -input.back();
        for (std::size_t j = 0; j + 1 < n; ++j) {
            const T angle = std::numbers::pi_v<T> * (static_cast<T>(j) + T{1}) *
                            (static_cast<T>(k) + T{0.5}) / static_cast<T>(n);
            sum += T{2} * input[j] * std::sin(angle);
        }
        output[k] = sum;
    }
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> type4(std::span<const T> input) {
    const std::size_t n = input.size();
    std::vector<T> output(n, T{});
    if (n == 0) return output;

    for (std::size_t k = 0; k < n; ++k) {
        T sum{};
        for (std::size_t j = 0; j < n; ++j) {
            const T angle = std::numbers::pi_v<T> * (static_cast<T>(j) + T{0.5}) *
                            (static_cast<T>(k) + T{0.5}) / static_cast<T>(n);
            sum += input[j] * std::sin(angle);
        }
        output[k] = T{2} * sum;
    }
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> type1_bluestein(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) return {};
    if (n == std::numeric_limits<std::size_t>::max()) {
        throw std::length_error("DST-I workspace size overflow");
    }
    const std::size_t period = detail::checked_twice(n + 1);
    std::vector<fft::Complex<T>> extension(period);
    for (std::size_t j = 0; j < n; ++j) {
        extension[j + 1] = {input[j], T{0}};
        extension[period - (j + 1)] = {-input[j], T{0}};
    }

    const auto spectrum = fft::bluestein<T>(extension);
    std::vector<T> output(n);
    for (std::size_t k = 0; k < n; ++k) output[k] = -spectrum[k + 1].imag();
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> type2_bluestein(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) return {};
    const std::size_t period = detail::checked_twice(n);
    std::vector<fft::Complex<T>> extension(period);
    for (std::size_t j = 0; j < n; ++j) {
        extension[j] = {input[j], T{0}};
        extension[period - 1 - j] = {-input[j], T{0}};
    }

    const auto spectrum = fft::bluestein<T>(extension);
    std::vector<T> output(n);
    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t m = k + 1;
        const T angle = std::numbers::pi_v<T> * static_cast<T>(m) /
                        (T{2} * static_cast<T>(n));
        const fft::Complex<T> phase{std::cos(angle), -std::sin(angle)};
        output[k] = -(spectrum[m] * phase).imag();
    }
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> type3_bluestein(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) return {};
    const std::size_t period = detail::checked_twice(n);
    std::vector<fft::Complex<T>> spectrum(period);

    for (std::size_t k = 0; k < n; ++k) {
        const std::size_t m = k + 1;
        if (m == n) {
            spectrum[m] = {input[k], T{0}};
            continue;
        }
        const T angle = std::numbers::pi_v<T> * static_cast<T>(m) /
                        (T{2} * static_cast<T>(n));
        const fft::Complex<T> phase{std::cos(angle), std::sin(angle)};
        const fft::Complex<T> value = fft::Complex<T>{T{0}, T{-1}} * input[k] * phase;
        spectrum[m] = value;
        spectrum[period - m] = std::conj(value);
    }

    const auto time = fft::bluestein<T>(spectrum, fft::Direction::inverse);
    std::vector<T> output(n);
    const T scale = static_cast<T>(period);
    for (std::size_t j = 0; j < n; ++j) output[j] = scale * time[j].real();
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> type4_bluestein(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) return {};
    const std::size_t extended_size = detail::checked_twice(n);
    std::vector<T> extended(extended_size, T{});
    for (std::size_t j = 0; j < n; ++j) extended[j] = input[j];
    const auto transformed = type2_bluestein<T>(extended);

    std::vector<T> output(n);
    for (std::size_t k = 0; k < n; ++k) output[k] = transformed[2 * k];
    return output;
}

}  // namespace signal_processing::dst
