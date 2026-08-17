#include "signal_processing/goertzel.hpp"
#include "signal_processing/hilbert.hpp"
#include "signal_processing/spectral.hpp"
#include "signal_processing/stft.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <numbers>
#include <vector>

namespace {

template <typename T>
[[nodiscard]] T tolerance() {
    if constexpr (std::same_as<T, float>) {
        return 2.0e-3F;
    } else {
        return 1.0e-10;
    }
}

template <typename T>
void expect_near(T actual, T expected, T multiplier = T{1}) {
    assert(std::abs(actual - expected) <=
           tolerance<T>() * multiplier * (T{1} + std::abs(expected)));
}

template <typename T>
void expect_near(std::complex<T> actual, std::complex<T> expected,
                 T multiplier = T{1}) {
    assert(std::abs(actual - expected) <=
           tolerance<T>() * multiplier * (T{1} + std::abs(expected)));
}

template <typename T>
void test_stft() {
    std::vector<T> input(19);
    for (std::size_t i = 0; i < input.size(); ++i) {
        input[i] = std::sin(static_cast<T>(0.31) * static_cast<T>(i)) +
                   static_cast<T>(0.1) * static_cast<T>(i);
    }
    std::vector<T> window(8, T{1});
    const auto transform = signal_processing::stft::bluestein<T>(input, window, 3, true);
    const auto restored = signal_processing::stft::inverse_bluestein<T>(transform, window);
    assert(restored.size() == input.size());
    for (std::size_t i = 0; i < input.size(); ++i)
        expect_near(restored[i], input[i], T{20});

    const auto power = signal_processing::stft::power_spectrogram<T>(transform, true);
    assert(power.size() == transform.frames.size());
    assert(power.front().size() == window.size() / 2 + 1);
}

template <typename T>
void test_complex_paths() {
    constexpr std::size_t n = 11;
    std::vector<std::complex<T>> input(n);
    for (std::size_t i = 0; i < n; ++i) {
        input[i] = {static_cast<T>(i % 4) - T{1},
                    static_cast<T>((i * 3) % 5) - T{2}};
    }

    std::vector<T> window(5, T{1});
    const auto transform =
        signal_processing::stft::bluestein<std::complex<T>>(input, window, 3, true);
    const auto restored =
        signal_processing::stft::inverse_bluestein<std::complex<T>>(transform, window);
    for (std::size_t i = 0; i < n; ++i) expect_near(restored[i], input[i], T{30});

    std::vector<T> full_window(n, T{1});
    const auto density =
        signal_processing::spectral::periodogram_bluestein<std::complex<T>>(
            input, full_window, static_cast<T>(n),
            signal_processing::spectral::Sides::two_sided);
    assert(density.size() == n);

    const auto coefficient = signal_processing::goertzel::bin<std::complex<T>>(input, 4);
    const auto spectrum = signal_processing::fft::bluestein<T>(input);
    expect_near(coefficient, spectrum[4], T{70});
}

template <typename T>
void test_periodogram_welch_and_csd() {
    constexpr std::size_t n = 16;
    constexpr std::size_t tone_bin = 3;
    std::vector<T> input(n);
    for (std::size_t i = 0; i < n; ++i) {
        input[i] = std::cos(T{2} * std::numbers::pi_v<T> *
                            static_cast<T>(tone_bin * i) / static_cast<T>(n));
    }
    std::vector<T> window(n, T{1});
    const auto periodogram = signal_processing::spectral::periodogram_bluestein<T>(
        input, window, static_cast<T>(n), signal_processing::spectral::Sides::one_sided);

    std::size_t peak = 0;
    for (std::size_t k = 1; k < periodogram.size(); ++k)
        if (periodogram[k] > periodogram[peak]) peak = k;
    assert(peak == tone_bin);

    T integrated{};
    for (const T value : periodogram) integrated += value;  // df = fs / N = 1.
    expect_near(integrated, T{0.5}, T{20});

    const auto welch = signal_processing::spectral::welch_bluestein<T>(
        input, window, n, static_cast<T>(n), signal_processing::spectral::Sides::one_sided);
    for (std::size_t k = 0; k < periodogram.size(); ++k)
        expect_near(welch[k], periodogram[k], T{10});

    const auto csd = signal_processing::spectral::cross_spectral_density_bluestein<T>(
        input, input, window, n, static_cast<T>(n),
        signal_processing::spectral::Sides::one_sided);
    for (std::size_t k = 0; k < periodogram.size(); ++k) {
        expect_near(csd[k].real(), periodogram[k], T{10});
        expect_near(csd[k].imag(), T{}, T{10});
    }

    const auto spectrogram = signal_processing::spectral::spectrogram_bluestein<T>(
        input, window, n, static_cast<T>(n), signal_processing::spectral::Sides::one_sided,
        false);
    assert(spectrogram.size() == 1);
    for (std::size_t k = 0; k < periodogram.size(); ++k)
        expect_near(spectrogram.front()[k], periodogram[k], T{10});

    std::vector<std::complex<T>> simple{{T{3}, T{4}}, {T{-1}, T{0}}};
    const auto magnitudes = signal_processing::spectral::magnitude<T>(simple);
    const auto powers = signal_processing::spectral::power<T>(simple);
    const auto phases = signal_processing::spectral::phase<T>(simple);
    expect_near(magnitudes[0], T{5});
    expect_near(powers[0], T{25});
    expect_near(phases[1], std::numbers::pi_v<T>, T{4});
}

template <typename T>
void test_goertzel() {
    constexpr std::size_t n = 15;
    std::vector<T> input(n);
    for (std::size_t i = 0; i < n; ++i)
        input[i] = static_cast<T>((i * 7) % 9) - T{4};

    std::vector<std::complex<T>> complex_input(n);
    for (std::size_t i = 0; i < n; ++i) complex_input[i] = {input[i], T{}};
    const auto spectrum = signal_processing::fft::bluestein<T>(complex_input);
    for (std::size_t k = 0; k < n; ++k)
        expect_near(signal_processing::goertzel::bin<T>(input, k), spectrum[k], T{50});

    const T frequency = static_cast<T>(0.137);
    std::complex<T> expected{};
    for (std::size_t i = 0; i < n; ++i) {
        const T angle = -T{2} * std::numbers::pi_v<T> * frequency * static_cast<T>(i);
        expected += input[i] * std::complex<T>{std::cos(angle), std::sin(angle)};
    }
    expect_near(signal_processing::goertzel::frequency<T>(input, frequency), expected, T{70});
}

template <typename T>
void test_hilbert() {
    constexpr std::size_t n = 32;
    constexpr std::size_t tone_bin = 5;
    std::vector<T> input(n);
    for (std::size_t i = 0; i < n; ++i) {
        input[i] = std::cos(T{2} * std::numbers::pi_v<T> *
                            static_cast<T>(tone_bin * i) / static_cast<T>(n));
    }

    const auto analytic = signal_processing::hilbert::analytic_bluestein<T>(input);
    const auto transformed = signal_processing::hilbert::transform_bluestein<T>(input);
    const auto envelope = signal_processing::hilbert::envelope_bluestein<T>(input);
    for (std::size_t i = 0; i < n; ++i) {
        const T phase = T{2} * std::numbers::pi_v<T> * static_cast<T>(tone_bin * i) /
                        static_cast<T>(n);
        expect_near(analytic[i].real(), std::cos(phase), T{30});
        expect_near(analytic[i].imag(), std::sin(phase), T{30});
        expect_near(transformed[i], std::sin(phase), T{30});
        expect_near(envelope[i], T{1}, T{30});
    }
}

}  // namespace

int main() {
    test_stft<float>();
    test_stft<double>();
    test_complex_paths<float>();
    test_complex_paths<double>();
    test_periodogram_welch_and_csd<float>();
    test_periodogram_welch_and_csd<double>();
    test_goertzel<float>();
    test_goertzel<double>();
    test_hilbert<float>();
    test_hilbert<double>();
}
