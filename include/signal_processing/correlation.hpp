#pragma once

#include "signal_processing/convolution.hpp"
#include "signal_processing/detail/sample.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_processing::correlation {

namespace detail {

struct Overlap {
    std::size_t lhs_first{};
    std::size_t rhs_first{};
    std::size_t count{};
};

[[nodiscard]] inline Overlap overlap_at(std::size_t lhs_size, std::size_t rhs_size,
                                        std::size_t output_index) {
    if (lhs_size == 0 || rhs_size == 0) return {};
    const std::size_t zero_lag_index = rhs_size - 1;
    if (output_index >= zero_lag_index) {
        const std::size_t lag = output_index - zero_lag_index;
        if (lag >= lhs_size) return {};
        return {lag, 0, std::min(lhs_size - lag, rhs_size)};
    }
    const std::size_t lag_magnitude = zero_lag_index - output_index;
    if (lag_magnitude >= rhs_size) return {};
    return {0, lag_magnitude, std::min(lhs_size, rhs_size - lag_magnitude)};
}

template <signal_processing::detail::Sample T>
[[nodiscard]] std::vector<signal_processing::detail::scalar_t<T>> energy_prefix(
    std::span<const T> input) {
    using Scalar = signal_processing::detail::scalar_t<T>;
    if (input.size() == std::numeric_limits<std::size_t>::max())
        throw std::length_error("correlation energy prefix size overflows size_t");
    std::vector<Scalar> prefix(input.size() + 1, Scalar{});
    for (std::size_t i = 0; i < input.size(); ++i)
        prefix[i + 1] = prefix[i] + signal_processing::detail::magnitude_squared(input[i]);
    return prefix;
}

template <fft::Scalar T>
[[nodiscard]] inline T range_sum(const std::vector<T>& prefix, std::size_t first,
                                 std::size_t count) {
    return prefix[first + count] - prefix[first];
}

template <signal_processing::detail::Sample T>
[[nodiscard]] std::vector<T> normalize(std::span<const T> lhs, std::span<const T> rhs,
                                       std::span<const T> raw) {
    using Scalar = signal_processing::detail::scalar_t<T>;
    const std::size_t expected = convolution::detail::linear_size(lhs.size(), rhs.size());
    if (raw.size() != expected) throw std::invalid_argument("correlation raw size mismatch");
    if (raw.empty()) return {};

    const auto lhs_prefix = energy_prefix<T>(lhs);
    const auto rhs_prefix = energy_prefix<T>(rhs);
    std::vector<T> output(raw.size(), T{});
    for (std::size_t index = 0; index < raw.size(); ++index) {
        const auto overlap = overlap_at(lhs.size(), rhs.size(), index);
        if (overlap.count == 0) continue;
        const Scalar lhs_energy = range_sum(lhs_prefix, overlap.lhs_first, overlap.count);
        const Scalar rhs_energy = range_sum(rhs_prefix, overlap.rhs_first, overlap.count);
        const Scalar denominator = std::sqrt(lhs_energy * rhs_energy);
        if (denominator != Scalar{}) output[index] = raw[index] / denominator;
    }
    return output;
}

}  // namespace detail

// r_xy[lag] = sum_n x[n] conj(y[n-lag]); output index = lag + rhs.size() - 1.
// https://en.wikipedia.org/wiki/Cross-correlation
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> cross_direct(std::span<const T> lhs,
                                                  std::span<const T> rhs) {
    const std::size_t size = convolution::detail::linear_size(lhs.size(), rhs.size());
    if (size == 0) return {};
    std::vector<T> output(size, T{});
    for (std::size_t index = 0; index < size; ++index) {
        const auto overlap = detail::overlap_at(lhs.size(), rhs.size(), index);
        T sum{};
        for (std::size_t offset = 0; offset < overlap.count; ++offset) {
            sum += lhs[overlap.lhs_first + offset] * signal_processing::detail::conjugate(
                       rhs[overlap.rhs_first + offset]);
        }
        output[index] = sum;
    }
    return output;
}

// Cross-correlation through convolution with the reversed conjugate of rhs.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> cross_fft(std::span<const T> lhs,
                                               std::span<const T> rhs) {
    if (lhs.empty() || rhs.empty()) return {};
    std::vector<T> reversed(rhs.size());
    for (std::size_t i = 0; i < rhs.size(); ++i)
        reversed[i] = signal_processing::detail::conjugate(rhs[rhs.size() - 1 - i]);
    return convolution::fft<T>(lhs, reversed);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> auto_direct(std::span<const T> input) {
    return cross_direct<T>(input, input);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> auto_fft(std::span<const T> input) {
    return cross_fft<T>(input, input);
}

// Per-lag normalized cross-correlation. Zero-energy overlaps map to zero.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> normalized_direct(std::span<const T> lhs,
                                                       std::span<const T> rhs) {
    const auto raw = cross_direct<T>(lhs, rhs);
    return detail::normalize<T>(lhs, rhs, raw);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> normalized_fft(std::span<const T> lhs,
                                                    std::span<const T> rhs) {
    const auto raw = cross_fft<T>(lhs, rhs);
    return detail::normalize<T>(lhs, rhs, raw);
}

// Compatibility names from the initial milestone; both are direct definitions.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> cross(std::span<const T> lhs, std::span<const T> rhs) {
    return cross_direct<T>(lhs, rhs);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> auto_correlation(std::span<const T> input) {
    return auto_direct<T>(input);
}

}  // namespace signal_processing::correlation
