#pragma once

#include "signal_processing/detail/sample.hpp"
#include "signal_processing/fft.hpp"
#include "signal_processing/stft.hpp"

#include <cmath>
#include <complex>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_processing::spectral {

enum class Sides { two_sided, one_sided };

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> magnitude(std::span<const fft::Complex<T>> spectrum) {
    std::vector<T> output(spectrum.size());
    for (std::size_t i = 0; i < spectrum.size(); ++i) output[i] = std::abs(spectrum[i]);
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> phase(std::span<const fft::Complex<T>> spectrum) {
    std::vector<T> output(spectrum.size());
    for (std::size_t i = 0; i < spectrum.size(); ++i) output[i] = std::arg(spectrum[i]);
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> power(std::span<const fft::Complex<T>> spectrum) {
    std::vector<T> output(spectrum.size());
    for (std::size_t i = 0; i < spectrum.size(); ++i) output[i] = std::norm(spectrum[i]);
    return output;
}

namespace detail {

template <signal_processing::detail::Sample Sample>
void validate_sides(Sides sides) {
    if (sides == Sides::one_sided && signal_processing::detail::SampleTraits<Sample>::is_complex)
        throw std::invalid_argument("one-sided spectral density is defined here only for real input");
}

template <fft::Scalar T>
[[nodiscard]] T window_energy(std::span<const T> window) {
    T energy{};
    for (const T value : window) energy += value * value;
    if (!(energy > T{})) throw std::invalid_argument("spectral window must have positive energy");
    return energy;
}

[[nodiscard]] inline std::size_t bins(std::size_t n, Sides sides) {
    return sides == Sides::one_sided ? n / 2 + 1 : n;
}

[[nodiscard]] inline bool doubled_one_sided_bin(std::size_t k, std::size_t n) noexcept {
    if (k == 0) return false;
    if ((n & 1U) == 0U && k == n / 2) return false;
    return k <= n / 2;
}

template <signal_processing::detail::Sample Sample>
[[nodiscard]] std::vector<signal_processing::detail::complex_t<Sample>> windowed_fft(
    std::span<const Sample> input,
    std::span<const signal_processing::detail::scalar_t<Sample>> window) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    using Complex = signal_processing::detail::complex_t<Sample>;
    if (input.size() != window.size())
        throw std::invalid_argument("spectral window size must match the input size");
    std::vector<Complex> time(input.size());
    for (std::size_t i = 0; i < input.size(); ++i)
        time[i] = signal_processing::detail::to_complex(input[i]) * window[i];
    return fft::bluestein<Scalar>(time, fft::Direction::forward);
}

template <fft::Scalar T>
void apply_one_sided_density_scaling(std::span<T> values, std::size_t transform_size) {
    for (std::size_t k = 0; k < values.size(); ++k)
        if (doubled_one_sided_bin(k, transform_size)) values[k] *= T{2};
}

template <fft::Scalar T>
void apply_one_sided_density_scaling(std::span<fft::Complex<T>> values,
                                     std::size_t transform_size) {
    for (std::size_t k = 0; k < values.size(); ++k)
        if (doubled_one_sided_bin(k, transform_size)) values[k] *= T{2};
}

}  // namespace detail

// Periodogram power spectral density.
// https://en.wikipedia.org/wiki/Periodogram
template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline std::vector<signal_processing::detail::scalar_t<Sample>>
periodogram_bluestein(
    std::span<const Sample> input,
    std::span<const signal_processing::detail::scalar_t<Sample>> window,
    signal_processing::detail::scalar_t<Sample> sample_rate,
    Sides sides = Sides::two_sided) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    detail::validate_sides<Sample>(sides);
    if (!(sample_rate > Scalar{})) throw std::invalid_argument("sample rate must be positive");
    if (input.empty()) return {};

    const Scalar energy = detail::window_energy<Scalar>(window);
    const auto spectrum = detail::windowed_fft<Sample>(input, window);
    const std::size_t count = detail::bins(input.size(), sides);
    const Scalar scale = Scalar{1} / (sample_rate * energy);
    std::vector<Scalar> output(count);
    for (std::size_t k = 0; k < count; ++k) output[k] = std::norm(spectrum[k]) * scale;
    if (sides == Sides::one_sided)
        detail::apply_one_sided_density_scaling<Scalar>(output, input.size());
    return output;
}

// Welch's method.
// https://en.wikipedia.org/wiki/Welch%27s_method
template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline std::vector<signal_processing::detail::scalar_t<Sample>> welch_bluestein(
    std::span<const Sample> input,
    std::span<const signal_processing::detail::scalar_t<Sample>> window,
    std::size_t hop_size,
    signal_processing::detail::scalar_t<Sample> sample_rate,
    Sides sides = Sides::two_sided) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    detail::validate_sides<Sample>(sides);
    if (window.empty()) throw std::invalid_argument("Welch window must be nonempty");
    if (hop_size == 0) throw std::invalid_argument("Welch hop size must be nonzero");
    if (!(sample_rate > Scalar{})) throw std::invalid_argument("sample rate must be positive");
    if (input.size() < window.size()) return {};

    const std::size_t frame_count = 1 + (input.size() - window.size()) / hop_size;
    const std::size_t count = detail::bins(window.size(), sides);
    std::vector<Scalar> output(count, Scalar{});
    for (std::size_t frame = 0; frame < frame_count; ++frame) {
        const auto segment = input.subspan(frame * hop_size, window.size());
        const auto estimate = periodogram_bluestein<Sample>(segment, window, sample_rate, sides);
        for (std::size_t k = 0; k < count; ++k) output[k] += estimate[k];
    }
    const Scalar average = Scalar{1} / static_cast<Scalar>(frame_count);
    for (auto& value : output) value *= average;
    return output;
}

template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline std::vector<signal_processing::detail::complex_t<Sample>>
cross_spectral_density_bluestein(
    std::span<const Sample> lhs, std::span<const Sample> rhs,
    std::span<const signal_processing::detail::scalar_t<Sample>> window,
    std::size_t hop_size,
    signal_processing::detail::scalar_t<Sample> sample_rate,
    Sides sides = Sides::two_sided) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    using Complex = signal_processing::detail::complex_t<Sample>;
    detail::validate_sides<Sample>(sides);
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("cross-spectral inputs must have equal length");
    if (window.empty()) throw std::invalid_argument("CSD window must be nonempty");
    if (hop_size == 0) throw std::invalid_argument("CSD hop size must be nonzero");
    if (!(sample_rate > Scalar{})) throw std::invalid_argument("sample rate must be positive");
    if (lhs.size() < window.size()) return {};

    const Scalar energy = detail::window_energy<Scalar>(window);
    const Scalar scale = Scalar{1} / (sample_rate * energy);
    const std::size_t frame_count = 1 + (lhs.size() - window.size()) / hop_size;
    const std::size_t count = detail::bins(window.size(), sides);
    std::vector<Complex> output(count, Complex{});

    for (std::size_t frame = 0; frame < frame_count; ++frame) {
        const std::size_t offset = frame * hop_size;
        const auto left = detail::windowed_fft<Sample>(lhs.subspan(offset, window.size()), window);
        const auto right = detail::windowed_fft<Sample>(rhs.subspan(offset, window.size()), window);
        for (std::size_t k = 0; k < count; ++k)
            output[k] += left[k] * std::conj(right[k]) * scale;
    }
    const Scalar average = Scalar{1} / static_cast<Scalar>(frame_count);
    for (auto& value : output) value *= average;
    if (sides == Sides::one_sided)
        detail::apply_one_sided_density_scaling<Scalar>(output, window.size());
    return output;
}

template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline std::vector<std::vector<signal_processing::detail::scalar_t<Sample>>>
spectrogram_bluestein(
    std::span<const Sample> input,
    std::span<const signal_processing::detail::scalar_t<Sample>> window,
    std::size_t hop_size,
    signal_processing::detail::scalar_t<Sample> sample_rate,
    Sides sides = Sides::two_sided,
    bool pad_end = false) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    detail::validate_sides<Sample>(sides);
    if (!(sample_rate > Scalar{})) throw std::invalid_argument("sample rate must be positive");
    const Scalar energy = detail::window_energy<Scalar>(window);
    const auto transform = stft::bluestein<Sample>(input, window, hop_size, pad_end);
    const std::size_t count = detail::bins(window.size(), sides);
    const Scalar scale = Scalar{1} / (sample_rate * energy);
    std::vector<std::vector<Scalar>> output;
    output.reserve(transform.frames.size());
    for (const auto& frame : transform.frames) {
        std::vector<Scalar> row(count);
        for (std::size_t k = 0; k < count; ++k) row[k] = std::norm(frame[k]) * scale;
        if (sides == Sides::one_sided)
            detail::apply_one_sided_density_scaling<Scalar>(row, window.size());
        output.push_back(std::move(row));
    }
    return output;
}

}  // namespace signal_processing::spectral
