#include "signal_processing/signal_processing.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <numbers>
#include <span>
#include <type_traits>
#include <vector>

namespace {

template <typename T>
[[nodiscard]] T tolerance() {
    if constexpr (std::same_as<T, float>) return 2.0e-3F;
    return 2.0e-10;
}

template <typename T>
void expect_near(std::span<const T> actual, std::span<const T> expected, double multiplier = 1.0) {
    assert(actual.size() == expected.size());
    using std::abs;
    const double eps = static_cast<double>(tolerance<typename signal_processing::detail::SampleTraits<T>::scalar_type>()) * multiplier;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        const double scale = 1.0 + static_cast<double>(abs(expected[i]));
        assert(static_cast<double>(abs(actual[i] - expected[i])) <= eps * scale);
    }
}

template <typename T>
[[nodiscard]] std::vector<T> make_signal(std::size_t n) {
    std::vector<T> output(n);
    for (std::size_t i = 0; i < n; ++i) {
        const int value = static_cast<int>((7 * i + 5) % 17) - 8;
        if constexpr (signal_processing::detail::SampleTraits<T>::is_complex) {
            using Scalar = signal_processing::detail::scalar_t<T>;
            output[i] = T{static_cast<Scalar>(value) / Scalar{5},
                          static_cast<Scalar>(static_cast<int>(i % 5) - 2) / Scalar{7}};
        } else {
            output[i] = static_cast<T>(value) / T{5};
        }
    }
    return output;
}

template <typename T>
void test_fir_execution() {
    const auto input = make_signal<T>(31);
    const std::vector<T> coefficients{T{0.125}, T{0.25}, T{0.375}, T{0.25}};
    const auto reference = signal_processing::fir::direct<T>(input, coefficients);
    expect_near<T>(signal_processing::fir::fft<T>(input, coefficients), reference, 12.0);
    expect_near<T>(signal_processing::fir::overlap_save<T>(input, coefficients, 16), reference, 12.0);

    signal_processing::fir::Filter<T> streaming(coefficients);
    expect_near<T>(streaming.process(input), reference);
    streaming.reset();
    expect_near<T>(streaming.process(input), reference);

    signal_processing::fir::OverlapSaveFilter<T> blocked(coefficients, 16);
    const std::size_t block_size = blocked.block_size();
    std::vector<T> block(block_size, T{}), produced(block_size, T{}), joined;
    for (std::size_t offset = 0; offset < input.size(); offset += block_size) {
        std::fill(block.begin(), block.end(), T{});
        const std::size_t count = std::min(block_size, input.size() - offset);
        std::copy_n(input.begin() + static_cast<std::ptrdiff_t>(offset), count, block.begin());
        blocked.process_block(block, produced);
        joined.insert(joined.end(), produced.begin(), produced.begin() + static_cast<std::ptrdiff_t>(count));
    }
    expect_near<T>(joined, reference, 12.0);
}

template <typename T>
void test_linear_phase_execution() {
    const auto input = make_signal<T>(23);
    const std::vector<T> symmetric{T{0.125}, T{0.25}, T{0.25}, T{0.25}, T{0.125}};
    assert(signal_processing::fir::symmetric_coefficients<T>(symmetric));
    expect_near<T>(signal_processing::fir::symmetric_direct<T>(input, symmetric),
                   signal_processing::fir::direct<T>(input, symmetric));

    const std::vector<T> antisymmetric{T{0.25}, T{0.5}, T{0}, T{-0.5}, T{-0.25}};
    assert(signal_processing::fir::antisymmetric_coefficients<T>(antisymmetric));
    expect_near<T>(signal_processing::fir::antisymmetric_direct<T>(input, antisymmetric),
                   signal_processing::fir::direct<T>(input, antisymmetric));
}

template <typename T>
[[nodiscard]] T response_magnitude(std::span<const T> coefficients, T frequency) {
    std::complex<T> response{};
    for (std::size_t n = 0; n < coefficients.size(); ++n) {
        const T angle = -T{2} * std::numbers::pi_v<T> * frequency * static_cast<T>(n);
        response += coefficients[n] * std::complex<T>{std::cos(angle), std::sin(angle)};
    }
    return std::abs(response);
}

template <typename T>
void test_fir_design() {
    const auto window = signal_processing::windows::hamming<T>(31);
    const auto lowpass = signal_processing::fir_design::lowpass_windowed_sinc<T>(31, T{0.15}, window);
    assert(signal_processing::fir::symmetric_coefficients<T>(lowpass, tolerance<T>()));
    assert(response_magnitude<T>(lowpass, T{0.05}) > response_magnitude<T>(lowpass, T{0.4}));

    const auto highpass = signal_processing::fir_design::highpass_windowed_sinc<T>(31, T{0.2}, window);
    assert(signal_processing::fir::symmetric_coefficients<T>(highpass, tolerance<T>()));
    assert(response_magnitude<T>(highpass, T{0.4}) > response_magnitude<T>(highpass, T{0.05}));

    const auto equiripple = signal_processing::fir_design::remez_lowpass<T>(31, T{0.18}, T{0.24}, T{1}, T{10});
    assert(signal_processing::fir::symmetric_coefficients<T>(equiripple, tolerance<T>()));
    T pass_error{};
    T stop_peak{};
    for (std::size_t i = 0; i <= 500; ++i) {
        const T f = T{0.5} * static_cast<T>(i) / T{500};
        const T magnitude = response_magnitude<T>(equiripple, f);
        if (f <= T{0.18}) pass_error = std::max(pass_error, std::abs(magnitude - T{1}));
        if (f >= T{0.24}) stop_peak = std::max(stop_peak, magnitude);
    }
    assert(pass_error > T{5} * stop_peak);
    assert(pass_error < T{15} * stop_peak);
}

template <typename T>
void test_multirate() {
    const auto input = make_signal<T>(13);
    const auto expanded = signal_processing::resampling::upsample_zero<T>(input, 3);
    assert(expanded.size() == (input.size() - 1) * 3 + 1);
    expect_near<T>(signal_processing::resampling::downsample<T>(expanded, 3), input);

    const auto held = signal_processing::resampling::zero_order_hold<T>(input, 3);
    for (std::size_t i = 0; i < input.size(); ++i)
        for (std::size_t p = 0; p < 3; ++p) assert(held[3 * i + p] == input[i]);

    const auto linear = signal_processing::resampling::linear<T>(input, 4);
    for (std::size_t i = 0; i < input.size(); ++i) assert(linear[4 * i] == input[i]);
    const auto sinc = signal_processing::resampling::windowed_sinc<T>(input, 4, 5);
    for (std::size_t i = 0; i < input.size(); ++i)
        assert(std::abs(sinc[4 * i] - input[i]) <= tolerance<signal_processing::detail::scalar_t<T>>() * 4);

    const std::vector<T> coefficients{T{0.05}, T{0.2}, T{0.5}, T{0.2}, T{0.05}};
    const auto direct = signal_processing::resampling::rational_fir<T>(input, 3, 2, coefficients);
    const auto polyphase = signal_processing::resampling::polyphase_rational<T>(input, 3, 2, coefficients);
    expect_near<T>(polyphase, direct, 4.0);

    const auto interpolation = signal_processing::resampling::interpolate_fir<T>(input, 3, coefficients);
    expect_near<T>(signal_processing::resampling::polyphase_interpolate<T>(input, 3, coefficients),
                   interpolation, 4.0);

    const auto decimation = signal_processing::resampling::decimate_fir<T>(input, 2, coefficients);
    expect_near<T>(signal_processing::resampling::polyphase_decimate<T>(input, 2, coefficients),
                   decimation, 4.0);
}

}  // namespace

int main() {
    test_fir_execution<float>();
    test_fir_execution<double>();
    test_fir_execution<std::complex<float>>();
    test_fir_execution<std::complex<double>>();
    test_linear_phase_execution<float>();
    test_linear_phase_execution<double>();
    test_fir_design<float>();
    test_fir_design<double>();
    test_multirate<float>();
    test_multirate<double>();
    test_multirate<std::complex<float>>();
    test_multirate<std::complex<double>>();
}
