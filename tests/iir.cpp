#include "signal_processing/signal_processing.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <span>
#include <type_traits>
#include <vector>

namespace iir = signal_processing::iir;
namespace design = signal_processing::iir::design;

namespace {

template <iir::Scalar T>
[[nodiscard]] T tolerance() {
    if constexpr (std::same_as<T, float>) {
        return 2.0e-3F;
    } else {
        return 2.0e-10;
    }
}

template <iir::Scalar T>
[[nodiscard]] std::vector<T> make_signal(std::size_t count) {
    std::vector<T> signal(count);
    for (std::size_t i = 0; i < count; ++i) {
        const T x = static_cast<T>(i);
        signal[i] = std::sin(static_cast<T>(0.17) * x) +
                    static_cast<T>(0.2) * std::cos(static_cast<T>(0.31) * x);
    }
    return signal;
}

template <iir::Scalar T>
void expect_near(T actual, T expected, T multiplier = T{1}) {
    assert(std::abs(actual - expected) <=
           tolerance<T>() * multiplier * (T{1} + std::abs(expected)));
}

template <iir::Scalar T>
void test_structures() {
    const std::vector<T> b{static_cast<T>(0.2), static_cast<T>(0.1), static_cast<T>(0.05)};
    const std::vector<T> a{T{1}, static_cast<T>(-0.4), static_cast<T>(0.12)};
    iir::DirectFormI<T> df1(b, a);
    iir::DirectFormII<T> df2(b, a);
    iir::TransposedDirectFormI<T> tdf1(b, a);
    iir::TransposedDirectFormII<T> tdf2(b, a);

    for (const T sample : make_signal<T>(128)) {
        const T reference = df1.process(sample);
        expect_near<T>(df2.process(sample), reference, T{4});
        expect_near<T>(tdf1.process(sample), reference, T{4});
        expect_near<T>(tdf2.process(sample), reference, T{4});
    }

    const iir::BiquadCoefficients<T> coefficients{
        static_cast<T>(0.2), static_cast<T>(0.1), static_cast<T>(0.05), T{1},
        static_cast<T>(-0.4), static_cast<T>(0.12)};
    iir::BiquadDirectFormI<T> biquad_df1(coefficients);
    iir::BiquadDirectFormII<T> biquad_df2(coefficients);
    iir::BiquadTransposedDirectFormII<T> biquad_tdf2(coefficients);
    df1.reset();
    for (const T sample : make_signal<T>(64)) {
        const T reference = df1.process(sample);
        expect_near<T>(biquad_df1.process(sample), reference, T{4});
        expect_near<T>(biquad_df2.process(sample), reference, T{4});
        expect_near<T>(biquad_tdf2.process(sample), reference, T{4});
    }
}

template <iir::Scalar T>
void test_prototypes_and_responses() {
    const T cutoff = static_cast<T>(0.25);
    const auto butterworth = design::butterworth<T>(4, design::Response::lowpass, cutoff);
    assert(design::stable(butterworth));
    expect_near<T>(std::abs(design::frequency_response(butterworth, cutoff)),
                   T{1} / std::sqrt(T{2}), T{12});

    const auto cheb1 = design::chebyshev1<T>(4, T{1}, design::Response::lowpass, cutoff);
    assert(design::stable(cheb1));
    expect_near<T>(std::abs(design::frequency_response(cheb1, cutoff)),
                   std::pow(T{10}, T{-1} / T{20}), T{16});

    const auto cheb2 = design::chebyshev2<T>(4, T{40}, design::Response::lowpass, cutoff);
    assert(design::stable(cheb2));
    expect_near<T>(std::abs(design::frequency_response(cheb2, cutoff)), static_cast<T>(0.01),
                   T{20});

    const auto elliptic =
        design::elliptic<T>(4, T{1}, T{40}, design::Response::lowpass, cutoff);
    assert(design::stable(elliptic));
    expect_near<T>(std::abs(design::frequency_response(elliptic, cutoff)),
                   std::pow(T{10}, T{-1} / T{20}), T{24});

    const auto bessel = design::bessel<T>(5, design::Response::lowpass, cutoff);
    assert(design::stable(bessel));

    assert(design::stable(
        design::butterworth<T>(3, design::Response::highpass, static_cast<T>(0.3))));
    assert(design::stable(design::butterworth<T>(
        3, design::Response::bandpass, static_cast<T>(0.2), static_cast<T>(0.4))));
    assert(design::stable(design::butterworth<T>(
        3, design::Response::bandstop, static_cast<T>(0.2), static_cast<T>(0.4))));
}

template <iir::Scalar T>
void test_sos_equivalence() {
    const auto zpk = design::elliptic<T>(5, T{1}, T{50}, design::Response::lowpass,
                                         static_cast<T>(0.28));
    const auto transfer = design::transfer_function(zpk);
    const auto sos = design::second_order_sections(zpk);
    iir::TransposedDirectFormII<T> direct(transfer.numerator, transfer.denominator);
    iir::SosCascade<T> cascade(sos.sections, sos.gain);

    for (std::size_t i = 0; i < 128; ++i) {
        const T sample = i == 0 ? T{1} : T{};
        expect_near<T>(cascade.process(sample), direct.process(sample), T{40});
    }
}

void test_analog_reference_points() {
    const auto elliptic = design::elliptic_prototype<double>(5, 1.0, 60.0);
    assert(elliptic.zeros.size() == 4);
    assert(elliptic.poles.size() == 5);
    assert(design::analog_stable(elliptic));
    expect_near<double>(elliptic.gain, 0.00750482523336344, 32.0);
    expect_near<double>(std::abs(design::analog_frequency_response(elliptic, 1.0)),
                        std::pow(10.0, -1.0 / 20.0), 32.0);

    const auto bessel =
        design::bessel_prototype<double>(6, design::BesselNormalization::magnitude_3db);
    assert(design::analog_stable(bessel));
    expect_near<double>(std::abs(design::analog_frequency_response(bessel, 1.0)),
                        1.0 / std::sqrt(2.0), 64.0);
}

void test_bessel_orders() {
    for (std::size_t order = 1; order <= 16; ++order) {
        const auto prototype =
            design::bessel_prototype<double>(order, design::BesselNormalization::delay);
        assert(prototype.poles.size() == order);
        for (const auto pole : prototype.poles) assert(pole.real() < 0.0);
    }
}

}  // namespace

int main() {
    static_assert(!iir::Scalar<long double>);
    test_structures<float>();
    test_structures<double>();
    test_prototypes_and_responses<float>();
    test_prototypes_and_responses<double>();
    test_sos_equivalence<float>();
    test_sos_equivalence<double>();
    test_bessel_orders();
    test_analog_reference_points();
}
