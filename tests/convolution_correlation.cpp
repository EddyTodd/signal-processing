#include "signal_processing/convolution.hpp"
#include "signal_processing/correlation.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace conv = signal_processing::convolution;
namespace corr = signal_processing::correlation;

template <typename T>
[[nodiscard]] T tolerance() {
    if constexpr (std::same_as<T, float>) return 4.0e-4F;
    else return 2.0e-10;
}

template <typename T>
void expect_near(
    std::span<const T> actual, std::span<const T> expected,
    signal_processing::detail::scalar_t<T> scale = signal_processing::detail::scalar_t<T>{1}) {
    using Scalar = signal_processing::detail::scalar_t<T>;
    assert(actual.size() == expected.size());
    const Scalar eps = tolerance<Scalar>() * scale;
    for (std::size_t i = 0; i < actual.size(); ++i)
        assert(std::abs(actual[i] - expected[i]) <= eps * (Scalar{1} + std::abs(expected[i])));
}

template <typename T>
std::vector<T> make_real_signal(std::size_t n) {
    std::vector<T> result(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto value = static_cast<int>((7 * i + 3) % 13) - 6;
        result[i] = static_cast<T>(value) / T{4};
    }
    return result;
}

template <typename T>
std::vector<std::complex<T>> make_complex_signal(std::size_t n) {
    std::vector<std::complex<T>> result(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto real = static_cast<int>((5 * i + 1) % 11) - 5;
        const auto imag = static_cast<int>((3 * i + 2) % 9) - 4;
        result[i] = {static_cast<T>(real) / T{3}, static_cast<T>(imag) / T{5}};
    }
    return result;
}

template <typename Sample>
void test_linear_catalog(const std::vector<Sample>& signal, const std::vector<Sample>& kernel) {
    using Scalar = signal_processing::detail::scalar_t<Sample>;
    const auto reference = conv::direct<Sample>(signal, kernel);
    expect_near<Sample>(conv::fft<Sample>(signal, kernel), reference, Scalar{8});
    expect_near<Sample>(conv::overlap_add<Sample>(signal, kernel, 3), reference, Scalar{10});

    std::size_t transform_size = 1;
    while (transform_size < kernel.size()) transform_size <<= 1;
    if (transform_size == kernel.size()) transform_size <<= 1;
    expect_near<Sample>(conv::overlap_save<Sample>(signal, kernel, transform_size), reference,
                        Scalar{10});
    expect_near<Sample>(conv::partitioned<Sample>(signal, kernel, 3), reference, Scalar{14});

    conv::StreamingDirect<Sample> streaming(kernel);
    auto streamed = streaming.process(signal);
    auto tail = streaming.flush();
    streamed.insert(streamed.end(), tail.begin(), tail.end());
    expect_near<Sample>(streamed, reference, Scalar{2});
}

template <typename T>
void test_real_catalog() {
    const auto signal = make_real_signal<T>(11);
    const auto kernel = make_real_signal<T>(5);
    test_linear_catalog(signal, kernel);

    const auto lhs = make_real_signal<T>(7);
    const auto rhs = make_real_signal<T>(7);
    expect_near<T>(conv::circular_bluestein<T>(lhs, rhs), conv::circular<T>(lhs, rhs), T{12});

    const auto direct_correlation = corr::cross_direct<T>(signal, kernel);
    expect_near<T>(corr::cross_fft<T>(signal, kernel), direct_correlation, T{10});
    expect_near<T>(corr::normalized_fft<T>(signal, kernel),
                   corr::normalized_direct<T>(signal, kernel), T{12});

    const auto normalized_auto = corr::normalized_direct<T>(signal, signal);
    assert(std::abs(normalized_auto[signal.size() - 1] - T{1}) <= tolerance<T>() * T{4});
}

template <typename T>
void test_complex_catalog() {
    using Complex = std::complex<T>;
    const auto signal = make_complex_signal<T>(9);
    const auto kernel = make_complex_signal<T>(4);
    test_linear_catalog(signal, kernel);

    const auto direct_correlation = corr::cross_direct<Complex>(signal, kernel);
    expect_near<Complex>(corr::cross_fft<Complex>(signal, kernel), direct_correlation, T{12});
    expect_near<Complex>(corr::normalized_fft<Complex>(signal, kernel),
                         corr::normalized_direct<Complex>(signal, kernel), T{14});
}

int main() {
    test_real_catalog<float>();
    test_real_catalog<double>();
    test_complex_catalog<float>();
    test_complex_catalog<double>();
}
