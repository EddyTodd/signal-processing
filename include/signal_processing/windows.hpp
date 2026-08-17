#pragma once

#include <cmath>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <vector>

namespace signal_processing::windows {

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> hann(std::size_t n) {
    std::vector<T> output(n, T{});
    if (n <= 1) {
        if (n == 1) output[0] = T{1};
        return output;
    }
    for (std::size_t i = 0; i < n; ++i) {
        output[i] = T{0.5} - T{0.5} * std::cos(T{2} * std::numbers::pi_v<T> *
                                               static_cast<T>(i) / static_cast<T>(n - 1));
    }
    return output;
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> hamming(std::size_t n) {
    std::vector<T> output(n, T{});
    if (n <= 1) {
        if (n == 1) output[0] = T{1};
        return output;
    }
    for (std::size_t i = 0; i < n; ++i) {
        output[i] = T{0.54} - T{0.46} * std::cos(T{2} * std::numbers::pi_v<T> *
                                                 static_cast<T>(i) / static_cast<T>(n - 1));
    }
    return output;
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> blackman(std::size_t n) {
    std::vector<T> output(n, T{});
    if (n <= 1) {
        if (n == 1) output[0] = T{1};
        return output;
    }
    for (std::size_t i = 0; i < n; ++i) {
        const T phase = T{2} * std::numbers::pi_v<T> * static_cast<T>(i) /
                        static_cast<T>(n - 1);
        output[i] = T{0.42} - T{0.5} * std::cos(phase) + T{0.08} * std::cos(T{2} * phase);
    }
    return output;
}

}  // namespace signal_processing::windows
