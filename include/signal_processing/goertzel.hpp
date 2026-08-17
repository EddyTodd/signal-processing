#pragma once

#include "signal_processing/detail/sample.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <stdexcept>

namespace signal_processing::goertzel {

// Goertzel algorithm.
// https://en.wikipedia.org/wiki/Goertzel_algorithm

template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline signal_processing::detail::complex_t<Sample> frequency(
    std::span<const Sample> input,
    signal_processing::detail::scalar_t<Sample> cycles_per_sample) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    using Complex = signal_processing::detail::complex_t<Sample>;

    if (!(cycles_per_sample >= Scalar{} && cycles_per_sample < Scalar{1}))
        throw std::invalid_argument("Goertzel frequency must be in [0,1) cycles/sample");
    if (input.empty()) return {};

    const Scalar omega = Scalar{2} * std::numbers::pi_v<Scalar> * cycles_per_sample;
    const Scalar cosine = std::cos(omega);
    const Scalar coefficient = Scalar{2} * cosine;
    Complex previous{};
    Complex previous2{};
    for (const auto& sample : input) {
        const Complex current = signal_processing::detail::to_complex(sample) +
                                coefficient * previous - previous2;
        previous2 = previous;
        previous = current;
    }

    // The recurrence state is referenced to the end of the block. Rotate it back
    // to the global DFT phase convention so the complex coefficient, not only its
    // magnitude, agrees with the forward DFT sum.
    const Scalar last_phase = -omega * static_cast<Scalar>(input.size() - 1);
    const Complex rotate_last{std::cos(last_phase), std::sin(last_phase)};
    const Scalar after_phase = last_phase - omega;
    const Complex rotate_after{std::cos(after_phase), std::sin(after_phase)};
    return rotate_last * previous - rotate_after * previous2;
}

template <signal_processing::detail::Sample Sample>
[[nodiscard]] inline signal_processing::detail::complex_t<Sample> bin(
    std::span<const Sample> input, std::size_t bin_index) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    if (input.empty()) {
        if (bin_index != 0) throw std::out_of_range("Goertzel bin is outside an empty transform");
        return {};
    }
    if (bin_index >= input.size()) throw std::out_of_range("Goertzel bin is outside the DFT");
    return frequency<Sample>(input, static_cast<Scalar>(bin_index) /
                                        static_cast<Scalar>(input.size()));
}

}  // namespace signal_processing::goertzel
