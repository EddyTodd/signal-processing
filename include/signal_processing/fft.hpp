#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_processing::fft {

template <typename T>
concept Scalar = std::floating_point<T>;

enum class Direction { forward, inverse };

template <Scalar T>
using Complex = std::complex<T>;

template <Scalar T>
[[nodiscard]] inline std::vector<Complex<T>> dft(std::span<const Complex<T>> input,
                                                 Direction direction = Direction::forward) {
    const std::size_t n = input.size();
    std::vector<Complex<T>> output(n);
    if (n == 0) {
        return output;
    }

    const T sign = direction == Direction::forward ? T{-1} : T{1};
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k = 0; k < n; ++k) {
        Complex<T> sum{};
        for (std::size_t t = 0; t < n; ++t) {
            const T angle = sign * tau * static_cast<T>(k) * static_cast<T>(t) /
                            static_cast<T>(n);
            sum += input[t] * Complex<T>{std::cos(angle), std::sin(angle)};
        }
        output[k] = sum;
    }

    if (direction == Direction::inverse) {
        const T scale = T{1} / static_cast<T>(n);
        for (auto& value : output) {
            value *= scale;
        }
    }
    return output;
}

[[nodiscard]] constexpr bool is_power_of_two(std::size_t n) noexcept {
    return n != 0 && (n & (n - 1)) == 0;
}

// Cooley-Tukey radix-2 FFT.
// https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
template <Scalar T>
inline void radix2_inplace(std::span<Complex<T>> data,
                           Direction direction = Direction::forward) {
    const std::size_t n = data.size();
    if (n <= 1) {
        return;
    }
    if (!is_power_of_two(n)) {
        throw std::invalid_argument("radix2_inplace requires a power-of-two length");
    }

    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(data[i], data[j]);
        }
    }

    const T sign = direction == Direction::forward ? T{-1} : T{1};
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t length = 2; length <= n;) {
        const T angle = sign * tau / static_cast<T>(length);
        const Complex<T> step{std::cos(angle), std::sin(angle)};
        for (std::size_t base = 0; base < n; base += length) {
            Complex<T> twiddle{T{1}, T{0}};
            for (std::size_t offset = 0; offset < length / 2; ++offset) {
                const Complex<T> upper = data[base + offset];
                const Complex<T> lower = data[base + offset + length / 2] * twiddle;
                data[base + offset] = upper + lower;
                data[base + offset + length / 2] = upper - lower;
                twiddle *= step;
            }
        }
        if (length == n) {
            break;
        }
        length <<= 1;
    }

    if (direction == Direction::inverse) {
        const T scale = T{1} / static_cast<T>(n);
        for (auto& value : data) {
            value *= scale;
        }
    }
}

template <Scalar T>
[[nodiscard]] inline std::vector<Complex<T>> radix2(
    std::span<const Complex<T>> input, Direction direction = Direction::forward) {
    std::vector<Complex<T>> output(input.begin(), input.end());
    radix2_inplace<T>(output, direction);
    return output;
}

}  // namespace signal_processing::fft
