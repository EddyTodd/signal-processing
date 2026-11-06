#pragma once

#include "signal_processing/detail/sample.hpp"
#include "signal_processing/fft.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace signal_processing::stft {

// Short-time Fourier transform.
// https://en.wikipedia.org/wiki/Short-time_Fourier_transform

template <signal_processing::detail::Sample Sample>
struct Result {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    using Complex = fft::Complex<Scalar>;

    std::size_t frame_size{};
    std::size_t hop_size{};
    std::size_t original_size{};
    bool padded_end{};
    std::vector<std::vector<Complex>> frames;
};

namespace detail {

// The forward transform accepts zero-valued windows (including an
// identically zero one) but cannot meaningfully process NaN/Inf weights.
// Reconstruction additionally needs every nonzero squared weight to
// survive accumulation in the public scalar format.
template <fft::Scalar T>
void validate_window(std::span<const T> window, bool reconstruction) {
    for (const T weight : window) {
        if (!std::isfinite(weight))
            throw std::invalid_argument("STFT window coefficients must be finite");
        if (reconstruction && weight != T{}) {
            const T squared = weight * weight;
            if (!(squared > T{}) || !std::isfinite(squared))
                throw std::invalid_argument(
                    "ISTFT nonzero window weights require finite nonzero squares");
        }
    }
}

[[nodiscard]] inline std::size_t frame_count(std::size_t signal_size, std::size_t frame_size,
                                             std::size_t hop_size, bool pad_end) {
    if (signal_size == 0) return 0;
    if (!pad_end) {
        if (signal_size < frame_size) return 0;
        return 1 + (signal_size - frame_size) / hop_size;
    }
    if (signal_size <= frame_size) return 1;
    const std::size_t remainder = signal_size - frame_size;
    const std::size_t required = 1 + remainder / hop_size +
                                 (remainder % hop_size == 0 ? 0 : 1);
    // When hop_size > frame_size, the regular hop grid can skip samples.
    // Never manufacture trailing frames whose start is past the input end:
    // such frames are entirely synthetic zeros, not signal-bearing windows.
    const std::size_t nonempty = 1 + (signal_size - 1) / hop_size;
    return std::min(required, nonempty);
}

}  // namespace detail

template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline Result<Sample> bluestein(
    std::span<const Sample> input,
    std::span<const signal_processing::detail::scalar_t<Sample>> window,
    std::size_t hop_size, bool pad_end = true) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    using Complex = fft::Complex<Scalar>;

    if (window.empty()) throw std::invalid_argument("STFT window must be nonempty");
    if (hop_size == 0) throw std::invalid_argument("STFT hop size must be nonzero");
    detail::validate_window<Scalar>(window, false);

    Result<Sample> result;
    result.frame_size = window.size();
    result.hop_size = hop_size;
    result.original_size = input.size();
    result.padded_end = pad_end;

    const std::size_t count = detail::frame_count(input.size(), window.size(), hop_size, pad_end);
    result.frames.reserve(count);
    std::vector<Complex> time(window.size());

    for (std::size_t frame = 0; frame < count; ++frame) {
        if (frame > std::numeric_limits<std::size_t>::max() / hop_size)
            throw std::length_error("STFT frame offset overflows size_t");
        const std::size_t offset = frame * hop_size;
        std::fill(time.begin(), time.end(), Complex{});
        const std::size_t available = offset < input.size() ? input.size() - offset : 0;
        const std::size_t samples = std::min(window.size(), available);
        for (std::size_t i = 0; i < samples; ++i) {
            time[i] = signal_processing::detail::to_complex(input[offset + i]) * window[i];
        }
        result.frames.push_back(fft::bluestein<Scalar>(time, fft::Direction::forward));
    }
    return result;
}

template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline std::vector<Sample> inverse_bluestein(
    const Result<Sample>& transform,
    std::span<const signal_processing::detail::scalar_t<Sample>> window,
    std::size_t output_size = 0) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    using Complex = fft::Complex<Scalar>;

    if (transform.frame_size == 0 || window.size() != transform.frame_size)
        throw std::invalid_argument("ISTFT window size must match the STFT frame size");
    if (transform.hop_size == 0) throw std::invalid_argument("ISTFT hop size must be nonzero");
    detail::validate_window<Scalar>(window, true);
    for (const auto& frame : transform.frames) {
        if (frame.size() != transform.frame_size)
            throw std::invalid_argument("ISTFT frame size mismatch");
        for (const auto bin : frame) {
            if (!std::isfinite(bin.real()) || !std::isfinite(bin.imag()))
                throw std::invalid_argument("ISTFT spectral frames must be finite");
        }
    }

    if (output_size == 0) output_size = transform.original_size;
    std::vector<Complex> accumulated(output_size, Complex{});
    std::vector<Scalar> normalization(output_size, Scalar{});

    for (std::size_t frame_index = 0; frame_index < transform.frames.size(); ++frame_index) {
        if (frame_index > std::numeric_limits<std::size_t>::max() / transform.hop_size)
            throw std::length_error("ISTFT frame offset overflows size_t");
        const std::size_t offset = frame_index * transform.hop_size;
        if (offset >= output_size) break;
        const auto time =
            fft::bluestein<Scalar>(transform.frames[frame_index], fft::Direction::inverse);
        const std::size_t count = std::min(transform.frame_size, output_size - offset);
        for (std::size_t i = 0; i < count; ++i) {
            accumulated[offset + i] += time[i] * window[i];
            normalization[offset + i] += window[i] * window[i];
        }
    }

    std::vector<Sample> output(output_size, Sample{});
    for (std::size_t i = 0; i < output_size; ++i) {
        if (normalization[i] != Scalar{}) accumulated[i] /= normalization[i];
        output[i] = signal_processing::detail::from_complex<Sample>(accumulated[i]);
    }
    return output;
}

template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline std::vector<std::vector<signal_processing::detail::scalar_t<Sample>>>
power_spectrogram(const Result<Sample>& transform, bool one_sided = false) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    if (one_sided && signal_processing::detail::SampleTraits<Sample>::is_complex)
        throw std::invalid_argument("one-sided STFT power is defined here only for real input");

    const std::size_t bins = one_sided ? transform.frame_size / 2 + 1 : transform.frame_size;
    std::vector<std::vector<Scalar>> output;
    output.reserve(transform.frames.size());
    for (const auto& frame : transform.frames) {
        if (frame.size() != transform.frame_size)
            throw std::invalid_argument("STFT frame size mismatch");
        std::vector<Scalar> row(bins);
        for (std::size_t k = 0; k < bins; ++k) row[k] = std::norm(frame[k]);
        output.push_back(std::move(row));
    }
    return output;
}

}  // namespace signal_processing::stft
