// Numerical contracts run with assertions ON even in Release.
#ifdef NDEBUG
#error "The numerical contract test must not be compiled with NDEBUG"
#endif
#include "signal_processing/fft.hpp"
#include "signal_processing/convolution.hpp"
#include "signal_processing/dct.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <span>
#include <utility>
#include <vector>

namespace {
namespace fft = signal_processing::fft;
namespace conv = signal_processing::convolution;
using C = std::complex<double>;
using Signal = std::vector<C>;

void close(C got, C expected, double tol = 2e-9) {
    assert(std::isfinite(got.real()) && std::isfinite(got.imag()));
    assert(std::abs(got - expected) <= tol * (1.0 + std::abs(expected)));
}
void close(double got, double expected, double tol = 2e-9) {
    assert(std::isfinite(got));
    assert(std::abs(got - expected) <= tol * (1.0 + std::abs(expected)));
}
void compare(std::span<const C> got, std::span<const C> expected, double tol = 2e-9) {
    assert(got.size() == expected.size());
    for (std::size_t i = 0; i < got.size(); ++i) close(got[i], expected[i], tol);
}
// Independent O(N^2) textbook DFT; no production FFT or DFT is called.
Signal reference_dft(std::span<const C> input) {
    const auto n = input.size();
    Signal output(n);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j) {
            const double angle = -2.0 * std::numbers::pi *
                static_cast<double>(j) * static_cast<double>(k) / static_cast<double>(n);
            output[k] += input[j] * C{std::cos(angle), std::sin(angle)};
        }
    return output;
}
Signal random_signal(std::size_t n, std::uint32_t seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    Signal result(n);
    for (auto& v : result) v = {dist(rng), dist(rng)};
    return result;
}
void transform_contracts(std::size_t n, std::uint32_t seed) {
    const auto input = random_signal(n, seed);
    const auto expected = reference_dft(input);
    compare(fft::dft<double>(input), expected);
    compare(fft::bluestein<double>(input), expected);
    compare(fft::mixed_radix<double>(input), expected);
    if (fft::is_power_of_two(n)) {
        compare(fft::radix2<double>(input), expected);
        compare(fft::stockham<double>(input), expected);
        compare(fft::split_radix<double>(input), expected);
    }
    const auto recovered = fft::bluestein<double>(expected, fft::Direction::inverse);
    compare(recovered, input);
    double input_energy = 0.0, output_energy = 0.0;
    for (const auto v : input) input_energy += std::norm(v);
    for (const auto v : expected) output_energy += std::norm(v);
    close(input_energy, output_energy / static_cast<double>(n), 5e-10);
    if (n > 1) {
        Signal shifted(n);
        for (std::size_t j = 0; j < n; ++j) shifted[(j + 1) % n] = input[j];
        const auto shifted_spectrum = fft::bluestein<double>(shifted);
        for (std::size_t k = 0; k < n; ++k) {
            const double phase = -2.0 * std::numbers::pi *
                static_cast<double>(k) / static_cast<double>(n);
            close(shifted_spectrum[k],
                  expected[k] * C{std::cos(phase),std::sin(phase)});
        }
    }
}
void convolution_contracts(std::size_t n, std::size_t m) {
    const auto a = random_signal(n, 811U + static_cast<std::uint32_t>(n));
    const auto b = random_signal(m, 419U + static_cast<std::uint32_t>(m));
    Signal reference(n + m - 1);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < m; ++j) reference[i+j] += a[i]*b[j];
    compare(conv::direct<C>(a,b), reference);
    compare(conv::fft<C>(a,b), reference, 5e-9);
    compare(conv::overlap_add<C>(a,b,7), reference, 5e-9);
}
void cosine_contracts(std::size_t n) {
    std::vector<double> samples(n);
    const auto seed = random_signal(n, 901U);
    for (std::size_t i = 0; i < n; ++i) samples[i] = seed[i].real();
    const auto expected = signal_processing::dct::type2<double>(samples);
    const auto got = signal_processing::dct::type2_bluestein<double>(samples);
    assert(expected.size() == got.size());
    for (std::size_t k = 0; k < n; ++k) close(got[k], expected[k]);
}
}  // namespace
int main() {
    for (const std::size_t n : {1U,2U,3U,5U,7U,8U,11U,16U,17U,31U,32U,
                                60U,64U,97U,127U,128U,255U,256U,257U})
        for (std::uint32_t seed : {0U,23U,4099U}) transform_contracts(n, seed);
    for (const auto [n,m] : std::array<std::pair<std::size_t,std::size_t>,7>{{
             {1,1},{2,17},{5,3},{17,9},{63,5},{128,17},{255,33}}})
        convolution_contracts(n,m);
    for (std::size_t n : {2U,3U,5U,16U,33U,127U}) cosine_contracts(n);
}
