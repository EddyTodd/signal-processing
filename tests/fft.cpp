#include "signal_processing/fft.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <type_traits>
#include <vector>

namespace fft = signal_processing::fft;

template <fft::Scalar T>
[[nodiscard]] T tolerance() {
    if constexpr (std::same_as<T, float>) {
        return 3.0e-4F;
    } else {
        return 2.0e-10;
    }
}

template <fft::Scalar T>
void expect_near(std::span<const std::complex<T>> actual,
                 std::span<const std::complex<T>> expected, T scale = T{1}) {
    assert(actual.size() == expected.size());
    const T eps = tolerance<T>() * scale;
    for (std::size_t i = 0; i < actual.size(); ++i)
        assert(std::abs(actual[i] - expected[i]) <= eps * (T{1} + std::abs(expected[i])));
}

template <fft::Scalar T>
std::vector<std::complex<T>> make_signal(std::size_t n) {
    std::vector<std::complex<T>> result(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto real = static_cast<int>((7 * i + 3) % 11) - 5;
        const auto imag = static_cast<int>((5 * i + 1) % 7) - 3;
        result[i] = {static_cast<T>(real) / T{3}, static_cast<T>(imag) / T{4}};
    }
    return result;
}

template <fft::Scalar T>
void test_power_of_two_catalog() {
    for (const std::size_t n : {1U, 2U, 4U, 8U, 16U, 32U}) {
        const auto input = make_signal<T>(n);
        const auto reference = fft::dft<T>(input);
        expect_near<T>(fft::radix2<T>(input), reference);
        expect_near<T>(fft::radix2_recursive<T>(input), reference, T{2});
        expect_near<T>(fft::stockham<T>(input), reference, T{2});
        expect_near<T>(fft::radix4<T>(input), reference, T{2});
        expect_near<T>(fft::split_radix<T>(input), reference, T{2});
        expect_near<T>(fft::modified_split_radix<T>(input), reference, T{12});

        const auto transformed = fft::split_radix<T>(input);
        expect_near<T>(fft::split_radix<T>(transformed, fft::Direction::inverse), input, T{3});
    }
}

template <fft::Scalar T>
void test_arbitrary_length_catalog() {
    for (const std::size_t n : {3U, 5U, 6U, 7U, 9U, 10U, 11U, 12U, 13U, 14U, 15U, 18U, 20U, 21U, 25U}) {
        const auto input = make_signal<T>(n);
        const auto reference = fft::dft<T>(input);
        expect_near<T>(fft::mixed_radix<T>(input), reference, T{6});
        expect_near<T>(fft::bluestein<T>(input), reference, T{20});
        expect_near<T>(fft::bluestein<T>(reference, fft::Direction::inverse), input, T{20});

        const auto factors = fft::detail::coprime_factor_split(n);
        if (factors.first > 1) {
            expect_near<T>(fft::good_thomas<T>(input, factors.first, factors.second), reference, T{10});
        }
        if (n == 3 || n == 5 || n == 7 || n == 11 || n == 13) {
            expect_near<T>(fft::rader<T>(input), reference, T{20});
            expect_near<T>(fft::rader<T>(reference, fft::Direction::inverse), input, T{20});
        }
    }
}

template <fft::Scalar T>
void test_codelets_and_plans() {
    for (const std::size_t radix : {2U, 3U, 4U, 5U, 7U}) {
        auto input = make_signal<T>(radix);
        const auto reference = fft::dft<T>(input);
        fft::SmallDftCodelet<T> codelet(radix);
        codelet.execute(input);
        expect_near<T>(input, reference, T{3});
    }

    auto x8 = make_signal<T>(8);
    const auto ref8 = fft::dft<T>(x8);
    fft::Radix2Plan<T> radix_plan(8);
    radix_plan.forward_inplace(x8);
    expect_near<T>(x8, ref8);
    radix_plan.inverse_inplace(x8);
    expect_near<T>(x8, make_signal<T>(8), T{2});

    auto x12 = make_signal<T>(12);
    const auto ref12 = fft::dft<T>(x12);
    fft::MixedRadixPlan<T> mixed_plan(12);
    std::vector<std::complex<T>> mixed_scratch(mixed_plan.scratch_size());
    mixed_plan.forward_inplace(x12, mixed_scratch);
    expect_near<T>(x12, ref12, T{6});
    mixed_plan.inverse_inplace(x12, mixed_scratch);
    expect_near<T>(x12, make_signal<T>(12), T{8});

    auto x15 = make_signal<T>(15);
    const auto ref15 = fft::dft<T>(x15);
    fft::GoodThomasPlan<T> pfa_plan(15);
    std::vector<std::complex<T>> pfa_scratch(pfa_plan.scratch_size());
    pfa_plan.forward_inplace(x15, pfa_scratch);
    expect_near<T>(x15, ref15, T{10});

    const auto x7 = make_signal<T>(7);
    const auto ref7 = fft::dft<T>(x7);
    fft::BluesteinPlan<T> bluestein_plan(7);
    std::vector<std::complex<T>> bluestein_scratch(bluestein_plan.scratch_size());
    std::vector<std::complex<T>> bluestein_output(7);
    bluestein_plan.forward(x7, bluestein_output, bluestein_scratch);
    expect_near<T>(bluestein_output, ref7, T{20});

    fft::RaderPlan<T> rader_plan(7);
    std::vector<std::complex<T>> rader_scratch(rader_plan.scratch_size());
    std::vector<std::complex<T>> rader_output(7);
    rader_plan.forward(x7, rader_output, rader_scratch);
    expect_near<T>(rader_output, ref7, T{20});
}

template <fft::Scalar T>
void test_real_plan() {
    constexpr std::size_t n = 16;
    std::vector<T> input(n);
    for (std::size_t i = 0; i < n; ++i) input[i] = static_cast<T>(static_cast<int>(i % 7) - 3) / T{2};
    fft::RealRadix2Plan<T> plan(n);
    std::vector<std::complex<T>> spectrum(plan.spectrum_size());
    std::vector<std::complex<T>> scratch(plan.scratch_size());
    plan.forward(input, spectrum, scratch);
    std::vector<T> restored(n);
    plan.inverse(spectrum, restored, scratch);
    const T eps = tolerance<T>() * T{4};
    for (std::size_t i = 0; i < n; ++i) assert(std::abs(restored[i] - input[i]) <= eps);
}

void test_explicit_kernels() {
    auto input = make_signal<double>(32);
    const auto reference = fft::dft<double>(input);

    fft::KernelRadix2Plan scalar(32, fft::KernelIsa::scalar);
    auto scalar_output = input;
    scalar.forward_inplace(scalar_output);
    expect_near<double>(scalar_output, reference);

    const auto capabilities = fft::kernel_capabilities();
    if (capabilities.avx2) {
        fft::KernelRadix2Plan avx2(32, fft::KernelIsa::avx2);
        auto output = input;
        avx2.forward_inplace(output);
        expect_near<double>(output, reference);
    }
    if (capabilities.avx512) {
        fft::KernelRadix2Plan avx512(32, fft::KernelIsa::avx512);
        auto output = input;
        avx512.forward_inplace(output);
        expect_near<double>(output, reference);
    }
}

int main() {
    static_assert(!fft::Scalar<long double>);
    test_power_of_two_catalog<float>();
    test_power_of_two_catalog<double>();
    test_arbitrary_length_catalog<float>();
    test_arbitrary_length_catalog<double>();
    test_codelets_and_plans<float>();
    test_codelets_and_plans<double>();
    test_real_plan<float>();
    test_real_plan<double>();
    test_explicit_kernels();
}
