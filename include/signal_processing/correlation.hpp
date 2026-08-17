#pragma once

#include "signal_processing/convolution.hpp"
#include "signal_processing/detail/sample.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

namespace signal_processing::correlation {

namespace detail {

template <signal_processing::detail::Sample T>
[[nodiscard]] std::vector<signal_processing::detail::scalar_t<T>> energy_prefix(
    std::span<const T> input) {
    using Scalar = signal_processing::detail::scalar_t<T>;
    std::vector<Scalar> prefix(input.size() + 1, Scalar{});
    for (std::size_t i = 0; i < input.size(); ++i)
        prefix[i + 1] = prefix[i] + signal_processing::detail::magnitude_squared(input[i]);
    return prefix;
}

template <fft::Scalar T>
[[nodiscard]] inline T range_sum(const std::vector<T>& prefix, std::size_t first,
                                 std::size_t last_exclusive) {
    return prefix[last_exclusive] - prefix[first];
}

template <signal_processing::detail::Sample T>
[[nodiscard]] std::vector<T> normalize(std::span<const T> lhs, std::span<const T> rhs,
                                       std::span<const T> raw) {
    using Scalar = signal_processing::detail::scalar_t<T>;
    const auto lhs_prefix = energy_prefix<T>(lhs);
    const auto rhs_prefix = energy_prefix<T>(rhs);
    std::vector<T> output(raw.size(), T{});
    const auto lhs_size = static_cast<std::ptrdiff_t>(lhs.size());
    const auto rhs_size = static_cast<std::ptrdiff_t>(rhs.size());

    for (std::size_t index = 0; index < raw.size(); ++index) {
        const auto lag = static_cast<std::ptrdiff_t>(index) - (rhs_size - 1);
        const auto i_first = std::max<std::ptrdiff_t>(0, lag);
        const auto i_last = std::min<std::ptrdiff_t>(lhs_size, rhs_size + lag);
        if (i_first >= i_last) continue;
        const auto j_first = i_first - lag;
        const auto j_last = i_last - lag;
        const Scalar lhs_energy = range_sum(lhs_prefix, static_cast<std::size_t>(i_first),
                                            static_cast<std::size_t>(i_last));
        const Scalar rhs_energy = range_sum(rhs_prefix, static_cast<std::size_t>(j_first),
                                            static_cast<std::size_t>(j_last));
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
    if (lhs.empty() || rhs.empty()) return {};
    const std::size_t size = convolution::detail::linear_size(lhs.size(), rhs.size());
    std::vector<T> output(size, T{});
    const auto rhs_extent = static_cast<std::ptrdiff_t>(rhs.size());
    for (std::size_t index = 0; index < size; ++index) {
        const auto lag = static_cast<std::ptrdiff_t>(index) - (rhs_extent - 1);
        T sum{};
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            const auto j = static_cast<std::ptrdiff_t>(i) - lag;
            if (j >= 0 && j < rhs_extent)
                sum += lhs[i] * signal_processing::detail::conjugate(
                                    rhs[static_cast<std::size_t>(j)]);
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
