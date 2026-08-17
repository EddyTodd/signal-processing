#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_processing::dct {

// Discrete cosine transform.
// https://en.wikipedia.org/wiki/Discrete_cosine_transform

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> type1(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) {
        return {};
    }
    if (n == 1) {
        return {input.front()};
    }
    std::vector<T> output(n, T{});
    for (std::size_t k = 0; k < n; ++k) {
        T sum = input.front() + ((k & 1U) == 0U ? input.back() : -input.back());
        for (std::size_t j = 1; j + 1 < n; ++j) {
            const T angle = std::numbers::pi_v<T> * static_cast<T>(j * k) /
                            static_cast<T>(n - 1);
            sum += T{2} * input[j] * std::cos(angle);
        }
        output[k] = sum;
    }
    return output;
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> type2(std::span<const T> input) {
    const std::size_t n = input.size();
    std::vector<T> output(n, T{});
    if (n == 0) {
        return output;
    }
    for (std::size_t k = 0; k < n; ++k) {
        T sum{};
        for (std::size_t j = 0; j < n; ++j) {
            const T angle = std::numbers::pi_v<T> *
                            (static_cast<T>(j) + T{0.5}) * static_cast<T>(k) /
                            static_cast<T>(n);
            sum += input[j] * std::cos(angle);
        }
        output[k] = T{2} * sum;
    }
    return output;
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> type3(std::span<const T> input) {
    const std::size_t n = input.size();
    std::vector<T> output(n, T{});
    if (n == 0) {
        return output;
    }
    for (std::size_t k = 0; k < n; ++k) {
        T sum = input.front();
        for (std::size_t j = 1; j < n; ++j) {
            const T angle = std::numbers::pi_v<T> * static_cast<T>(j) *
                            (static_cast<T>(k) + T{0.5}) / static_cast<T>(n);
            sum += T{2} * input[j] * std::cos(angle);
        }
        output[k] = sum;
    }
    return output;
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> type4(std::span<const T> input) {
    const std::size_t n = input.size();
    std::vector<T> output(n, T{});
    if (n == 0) {
        return output;
    }
    for (std::size_t k = 0; k < n; ++k) {
        T sum{};
        for (std::size_t j = 0; j < n; ++j) {
            const T angle = std::numbers::pi_v<T> *
                            (static_cast<T>(j) + T{0.5}) *
                            (static_cast<T>(k) + T{0.5}) / static_cast<T>(n);
            sum += input[j] * std::cos(angle);
        }
        output[k] = T{2} * sum;
    }
    return output;
}

}  // namespace signal_processing::dct
