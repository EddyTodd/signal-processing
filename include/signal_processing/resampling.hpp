#pragma once

#include "signal_processing/convolution.hpp"
#include "signal_processing/detail/sample.hpp"
#include "signal_processing/fir.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_processing::resampling {

// Downsampling / decimation / interpolation / rational sample-rate conversion.
// https://en.wikipedia.org/wiki/Downsampling_(signal_processing)
// https://en.wikipedia.org/wiki/Upsampling
// https://en.wikipedia.org/wiki/Polyphase_quadrature_filter

namespace detail {

[[nodiscard]] inline std::size_t checked_upsampled_size(std::size_t input_size,
                                                        std::size_t factor) {
    if (factor == 0) throw std::invalid_argument("resampling factor must be nonzero");
    if (input_size == 0) return 0;
    if (input_size - 1 > (std::numeric_limits<std::size_t>::max() - 1) / factor)
        throw std::length_error("upsampled size overflows size_t");
    return (input_size - 1) * factor + 1;
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline signal_processing::detail::scalar_t<T> sinc(
    signal_processing::detail::scalar_t<T> x) {
    using Scalar = signal_processing::detail::scalar_t<T>;
    if (x == Scalar{0}) return Scalar{1};
    const Scalar pix = std::numbers::pi_v<Scalar> * x;
    return std::sin(pix) / pix;
}

}  // namespace detail

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> downsample(std::span<const T> input, std::size_t factor,
                                                std::size_t phase = 0) {
    if (factor == 0) throw std::invalid_argument("downsample factor must be nonzero");
    if (phase >= factor) throw std::invalid_argument("downsample phase must be less than factor");
    if (phase >= input.size()) return {};
    const std::size_t count = 1 + (input.size() - 1 - phase) / factor;
    std::vector<T> output(count);
    for (std::size_t i = 0; i < count; ++i) output[i] = input[phase + i * factor];
    return output;
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> upsample_zero(std::span<const T> input,
                                                   std::size_t factor) {
    const std::size_t output_size = detail::checked_upsampled_size(input.size(), factor);
    std::vector<T> output(output_size, T{});
    for (std::size_t i = 0; i < input.size(); ++i) output[i * factor] = input[i];
    return output;
}

// Zero-order hold interpolation: repeat each input sample factor times.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> zero_order_hold(std::span<const T> input,
                                                     std::size_t factor) {
    if (factor == 0) throw std::invalid_argument("interpolation factor must be nonzero");
    if (input.size() > std::numeric_limits<std::size_t>::max() / factor)
        throw std::length_error("zero-order-hold size overflows size_t");
    std::vector<T> output(input.size() * factor);
    for (std::size_t i = 0; i < input.size(); ++i)
        std::fill_n(output.begin() + static_cast<std::ptrdiff_t>(i * factor), factor, input[i]);
    return output;
}

// Piecewise-linear interpolation with factor-1 new points between adjacent samples.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> linear(std::span<const T> input, std::size_t factor) {
    const std::size_t output_size = detail::checked_upsampled_size(input.size(), factor);
    if (output_size == 0) return {};
    using Scalar = signal_processing::detail::scalar_t<T>;
    std::vector<T> output(output_size);
    for (std::size_t segment = 0; segment + 1 < input.size(); ++segment) {
        for (std::size_t p = 0; p < factor; ++p) {
            const Scalar fraction = static_cast<Scalar>(p) / static_cast<Scalar>(factor);
            output[segment * factor + p] = input[segment] * (Scalar{1} - fraction) +
                                           input[segment + 1] * fraction;
        }
    }
    output.back() = input.back();
    return output;
}

// Finite Lanczos-windowed sinc interpolation. half_width is the sinc radius in input samples.
// https://en.wikipedia.org/wiki/Whittaker%E2%80%93Shannon_interpolation_formula
// https://en.wikipedia.org/wiki/Lanczos_resampling

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> windowed_sinc(std::span<const T> input, std::size_t factor,
                                                   std::size_t half_width = 4) {
    const std::size_t output_size = detail::checked_upsampled_size(input.size(), factor);
    if (output_size == 0) return {};
    if (half_width == 0) throw std::invalid_argument("sinc interpolation half-width must be nonzero");
    using Scalar = signal_processing::detail::scalar_t<T>;
    std::vector<T> output(output_size, T{});
    for (std::size_t m = 0; m < output_size; ++m) {
        const Scalar position = static_cast<Scalar>(m) / static_cast<Scalar>(factor);
        const auto center = static_cast<std::ptrdiff_t>(m / factor);
        const auto radius = static_cast<std::ptrdiff_t>(half_width);
        T sum{};
        for (std::ptrdiff_t n = center - radius + 1; n <= center + radius; ++n) {
            if (n < 0 || n >= static_cast<std::ptrdiff_t>(input.size())) continue;
            const Scalar distance = position - static_cast<Scalar>(n);
            if (std::abs(distance) >= static_cast<Scalar>(half_width)) continue;
            const Scalar weight = detail::sinc<T>(distance) *
                                  detail::sinc<T>(distance / static_cast<Scalar>(half_width));
            sum += input[static_cast<std::size_t>(n)] * weight;
        }
        output[m] = sum;
    }
    return output;
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> decimate_fir(std::span<const T> input,
                                                  std::size_t factor,
                                                  std::span<const T> coefficients,
                                                  std::size_t phase = 0) {
    return downsample<T>(fir::direct<T>(input, coefficients), factor, phase);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> interpolate_fir(std::span<const T> input,
                                                     std::size_t factor,
                                                     std::span<const T> coefficients) {
    const auto expanded = upsample_zero<T>(input, factor);
    return convolution::direct<T>(expanded, coefficients);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> rational_fir(std::span<const T> input,
                                                  std::size_t up_factor,
                                                  std::size_t down_factor,
                                                  std::span<const T> coefficients) {
    if (down_factor == 0) throw std::invalid_argument("downsample factor must be nonzero");
    const auto expanded = upsample_zero<T>(input, up_factor);
    const auto filtered = convolution::direct<T>(expanded, coefficients);
    return downsample<T>(filtered, down_factor);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<std::vector<T>> polyphase_decompose(
    std::span<const T> coefficients, std::size_t phase_count) {
    if (phase_count == 0) throw std::invalid_argument("polyphase count must be nonzero");
    std::vector<std::vector<T>> phases(phase_count);
    for (std::size_t k = 0; k < coefficients.size(); ++k)
        phases[k % phase_count].push_back(coefficients[k]);
    return phases;
}

// Rational L/M conversion evaluated directly in polyphase form. This is algebraically
// equivalent to zero-stuff by L, FIR filter, then keep every M-th output, but it never
// materializes or multiplies by the inserted zeros.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> polyphase_rational(std::span<const T> input,
                                                        std::size_t up_factor,
                                                        std::size_t down_factor,
                                                        std::span<const T> coefficients) {
    if (up_factor == 0 || down_factor == 0)
        throw std::invalid_argument("rational resampling factors must be nonzero");
    if (input.empty() || coefficients.empty()) return {};
    const std::size_t expanded_size = detail::checked_upsampled_size(input.size(), up_factor);
    if (expanded_size > std::numeric_limits<std::size_t>::max() - coefficients.size() + 1)
        throw std::length_error("rational FIR output size overflows size_t");
    const std::size_t full_size = expanded_size + coefficients.size() - 1;
    const std::size_t output_size = 1 + (full_size - 1) / down_factor;
    const auto phases = polyphase_decompose<T>(coefficients, up_factor);

    std::vector<T> output(output_size, T{});
    for (std::size_t r = 0; r < output_size; ++r) {
        if (r > std::numeric_limits<std::size_t>::max() / down_factor)
            throw std::length_error("rational resampling index overflows size_t");
        const std::size_t high_rate_index = r * down_factor;
        const std::size_t phase = high_rate_index % up_factor;
        const std::size_t input_center = high_rate_index / up_factor;
        const auto& phase_coefficients = phases[phase];
        T sum{};
        for (std::size_t q = 0; q < phase_coefficients.size(); ++q) {
            if (q > input_center) break;
            const std::size_t input_index = input_center - q;
            if (input_index < input.size()) sum += phase_coefficients[q] * input[input_index];
        }
        output[r] = sum;
    }
    return output;
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> polyphase_interpolate(std::span<const T> input,
                                                           std::size_t factor,
                                                           std::span<const T> coefficients) {
    return polyphase_rational<T>(input, factor, 1, coefficients);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> polyphase_decimate(std::span<const T> input,
                                                        std::size_t factor,
                                                        std::span<const T> coefficients) {
    return polyphase_rational<T>(input, 1, factor, coefficients);
}

}  // namespace signal_processing::resampling
