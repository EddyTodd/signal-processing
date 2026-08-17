#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>
#include <vector>

namespace signal_processing::correlation {

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> cross(std::span<const T> lhs, std::span<const T> rhs) {
    if (lhs.empty() || rhs.empty()) {
        return {};
    }
    const std::size_t size = lhs.size() + rhs.size() - 1;
    std::vector<T> output(size, T{});
    const auto rhs_extent = static_cast<std::ptrdiff_t>(rhs.size());
    for (std::size_t index = 0; index < size; ++index) {
        const auto lag = static_cast<std::ptrdiff_t>(index) - (rhs_extent - 1);
        T sum{};
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            const auto j = static_cast<std::ptrdiff_t>(i) - lag;
            if (j >= 0 && j < rhs_extent) {
                sum += lhs[i] * rhs[static_cast<std::size_t>(j)];
            }
        }
        output[index] = sum;
    }
    return output;
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> auto_correlation(std::span<const T> input) {
    return cross<T>(input, input);
}

}  // namespace signal_processing::correlation
