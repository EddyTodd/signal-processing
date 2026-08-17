#include "signal_processing/signal_processing.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <vector>

namespace {

bool near(double lhs, double rhs, double eps = 1e-9) {
    return std::abs(lhs - rhs) <= eps;
}

}  // namespace

int main() {
    using Complex = std::complex<double>;

    const std::vector<Complex> input{{1.0, 0.0}, {2.0, 0.0}, {3.0, 0.0}, {4.0, 0.0}};
    const auto direct_fft = signal_processing::fft::dft<double>(input);
    const auto radix2_fft = signal_processing::fft::radix2<double>(input);
    for (std::size_t i = 0; i < input.size(); ++i) {
        assert(std::abs(direct_fft[i] - radix2_fft[i]) < 1e-9);
    }
    const auto restored = signal_processing::fft::radix2<double>(
        radix2_fft, signal_processing::fft::Direction::inverse);
    for (std::size_t i = 0; i < input.size(); ++i) {
        assert(std::abs(restored[i] - input[i]) < 1e-9);
    }

    const std::vector<double> a{1.0, 2.0, 3.0};
    const std::vector<double> b{4.0, 5.0};
    const auto direct_conv = signal_processing::convolution::direct<double>(a, b);
    const auto fft_conv = signal_processing::convolution::fft<double>(a, b);
    assert(direct_conv.size() == 4);
    assert(near(direct_conv[0], 4.0));
    assert(near(direct_conv[1], 13.0));
    assert(near(direct_conv[2], 22.0));
    assert(near(direct_conv[3], 15.0));
    for (std::size_t i = 0; i < direct_conv.size(); ++i) {
        assert(near(direct_conv[i], fft_conv[i]));
    }

    const auto dct2 = signal_processing::dct::type2<double>(a);
    const auto dct3 = signal_processing::dct::type3<double>(dct2);
    for (std::size_t i = 0; i < a.size(); ++i) {
        assert(near(dct3[i] / (2.0 * static_cast<double>(a.size())), a[i]));
    }

    signal_processing::fir::Filter<double> fir({0.5, 0.5});
    assert(near(fir.process(2.0), 1.0));
    assert(near(fir.process(4.0), 3.0));

    signal_processing::iir::BiquadTransposedDirectFormII<double> identity(
        1.0, 0.0, 0.0, 1.0, 0.0, 0.0);
    assert(near(identity.process(3.5), 3.5));

    const auto autocorrelation = signal_processing::correlation::auto_correlation<double>(a);
    assert(autocorrelation.size() == 5);
    assert(near(autocorrelation[2], 14.0));

    const auto hann = signal_processing::windows::hann<double>(5);
    assert(near(hann.front(), 0.0));
    assert(near(hann[2], 1.0));
    assert(near(hann.back(), 0.0));
}
