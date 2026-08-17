#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_processing::mdct {

// Modified discrete cosine transform.
// https://en.wikipedia.org/wiki/Modified_discrete_cosine_transform

template <typename T>
concept Scalar = std::same_as<T, float> || std::same_as<T, double>;

template <Scalar T>
[[nodiscard]] inline std::vector<T> direct(std::span<const T> input) {
    if (input.empty()) return {};
    if ((input.size() & 1U) != 0U) {
        throw std::invalid_argument("MDCT requires an even input length 2N");
    }

    const std::size_t n = input.size() / 2;
    std::vector<T> output(n, T{});
    for (std::size_t k = 0; k < n; ++k) {
        T sum{};
        for (std::size_t j = 0; j < 2 * n; ++j) {
            const T angle = std::numbers::pi_v<T> / static_cast<T>(n) *
                            (static_cast<T>(j) + T{0.5} + static_cast<T>(n) / T{2}) *
                            (static_cast<T>(k) + T{0.5});
            sum += input[j] * std::cos(angle);
        }
        output[k] = sum;
    }
    return output;
}

template <Scalar T>
[[nodiscard]] inline std::vector<T> inverse_direct(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) return {};

    std::vector<T> output(2 * n, T{});
    const T scale = T{2} / static_cast<T>(n);
    for (std::size_t j = 0; j < 2 * n; ++j) {
        T sum{};
        for (std::size_t k = 0; k < n; ++k) {
            const T angle = std::numbers::pi_v<T> / static_cast<T>(n) *
                            (static_cast<T>(j) + T{0.5} + static_cast<T>(n) / T{2}) *
                            (static_cast<T>(k) + T{0.5});
            sum += input[k] * std::cos(angle);
        }
        output[j] = scale * sum;
    }
    return output;
}

}  // namespace signal_processing::mdct
