#include "signal_processing/signal_processing.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

template <typename Exception, typename Function>
[[nodiscard]] bool throws(Function&& function) {
    try {
        function();
    } catch (const Exception&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

int main() {
    namespace sp = signal_processing;

    static_assert(!sp::fft::Scalar<long double>);
    static_assert(!sp::windows::Scalar<long double>);
    static_assert(!sp::iir::Scalar<long double>);
    static_assert(!sp::detail::Sample<long double>);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    const std::vector<double> empty;
    const std::vector<double> signal{1.0, -2.0, 3.0};

    assert(throws<std::out_of_range>([&] { (void)sp::goertzel::bin<double>(empty, 0); }));

    assert(throws<std::invalid_argument>([&] {
        (void)sp::windows::kaiser<double>(8, infinity);
    }));
    assert(throws<std::invalid_argument>([&] {
        (void)sp::windows::gaussian<double>(8, nan);
    }));
    assert(throws<std::invalid_argument>([&] {
        (void)sp::windows::hann<double>(8, static_cast<sp::windows::Sampling>(99));
    }));

    assert(throws<std::invalid_argument>([&] { (void)sp::fir::direct<double>(signal, empty); }));
    assert(throws<std::invalid_argument>([&] {
        (void)sp::fir::fft<double>(std::span<const double>{}, empty);
    }));
    assert(!sp::fir::symmetric_coefficients<double>(empty));
    assert(throws<std::invalid_argument>([&] {
        (void)sp::fir::symmetric_coefficients<double>(signal, nan);
    }));
    assert(throws<std::invalid_argument>([&] {
        (void)sp::resampling::rational_fir<double>(signal, 2, 3, empty);
    }));

    assert(throws<std::length_error>([] {
        (void)sp::convolution::detail::linear_size(
            std::numeric_limits<std::size_t>::max(), 2);
    }));
    assert(throws<std::length_error>([] {
        (void)sp::fft::detail::next_power_of_two(std::numeric_limits<std::size_t>::max());
    }));
    assert(throws<std::invalid_argument>([] {
        (void)sp::fft::GoodThomasPlan<double>(6, 2, 4);
    }));

    const std::vector<std::complex<double>> fft_input{{1.0, 0.0}, {2.0, -1.0},
                                                       {-0.5, 2.0}, {3.0, 0.25},
                                                       {1.5, -0.5}, {-2.0, 1.0}};
    const auto pfa = sp::fft::good_thomas<double>(fft_input, sp::fft::Direction::forward, 2, 3);
    const auto reference = sp::fft::dft<double>(fft_input);
    assert(pfa.size() == reference.size());
    for (std::size_t i = 0; i < pfa.size(); ++i)
        assert(std::abs(pfa[i] - reference[i]) <= 1.0e-10 * (1.0 + std::abs(reference[i])));

    sp::stft::Result<double> malformed;
    malformed.frames.emplace_back();
    assert(throws<std::invalid_argument>([&] { (void)sp::stft::power_spectrogram(malformed); }));

    assert(throws<std::invalid_argument>([&] {
        (void)sp::iir::DirectFormI<double>({1.0}, {nan});
    }));
    assert(throws<std::invalid_argument>([&] {
        (void)sp::iir::BiquadTransposedDirectFormII<double>(1.0, 0.0, 0.0,
                                                            1.0, infinity, 0.0);
    }));
    sp::iir::DirectFormI<double> ieee_filter({1.0}, {1.0});
    assert(std::isnan(ieee_filter.process(nan)));

    assert(throws<std::invalid_argument>([] {
        (void)sp::iir::design::butterworth<double>(
            2, static_cast<sp::iir::design::Response>(99), 0.2);
    }));
    assert(throws<std::invalid_argument>([] {
        (void)sp::iir::design::bessel_prototype<double>(
            1, static_cast<sp::iir::design::BesselNormalization>(99));
    }));
    assert(throws<std::invalid_argument>([&] {
        (void)sp::iir::design::chebyshev1_prototype<double>(4, nan);
    }));
    const auto prototype = sp::iir::design::butterworth_prototype<double>(2);
    assert(throws<std::invalid_argument>([&] {
        (void)sp::iir::design::analog_frequency_response(prototype, nan);
    }));

    const std::vector<double> bands{0.0, 0.2, 0.3, 0.5};
    const std::vector<double> desired{nan, 0.0};
    const std::vector<double> weights{1.0, 1.0};
    assert(throws<std::invalid_argument>([&] {
        (void)sp::fir_design::remez_type1<double>(5, bands, desired, weights);
    }));

    return 0;
}
