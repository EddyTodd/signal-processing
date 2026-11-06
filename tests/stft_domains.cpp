// STFT parameter-domain and inverse reconstruction regression.
#ifdef NDEBUG
#error "STFT domains require assertions active in Release"
#endif

#include "signal_processing/stft.hpp"

#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
namespace stft = signal_processing::stft;

template <class Call>
void rejects_invalid(Call&& call) {
    bool rejected = false;
    try { (void)call(); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
}

template <class Sample>
void verify() {
    using T = signal_processing::detail::scalar_t<Sample>;
    const T nan = std::numeric_limits<T>::quiet_NaN();
    const T inf = std::numeric_limits<T>::infinity();
    const T tiny = std::numeric_limits<T>::denorm_min();
    const T huge = std::numeric_limits<T>::max() / T{2};
    const T tolerance = sizeof(T) == sizeof(float) ? static_cast<T>(0.002)
                                                   : static_cast<T>(3.e-10);

    std::vector<Sample> input(11);
    for (std::size_t i = 0; i < input.size(); ++i)
        input[i] = Sample{static_cast<T>(i % 7) - T{3}};
    const std::vector<T> window{T{1}, T{-1}, T{0.5}, T{0.75}};
    const auto transform = stft::bluestein<Sample>(input, window, 2, true);
    const auto restored = stft::inverse_bluestein<Sample>(transform, window);
    assert(restored.size() == input.size());
    for (std::size_t i = 0; i < input.size(); ++i) {
        assert(std::abs(restored[i] - input[i]) <=
               tolerance * (T{1} + std::abs(input[i])));
    }

    // Even if there are no signal-bearing frames, a window containing
    // NaN/Inf has an invalid transform parameter contract.
    const std::vector<Sample> empty;
    for (const T bad : {nan, inf, -inf}) {
        auto invalid = window;
        invalid[2] = bad;
        rejects_invalid([&] { return stft::bluestein<Sample>(
                                   empty, invalid, 2, true); });
        rejects_invalid([&] { return stft::bluestein<Sample>(
                                   input, invalid, 2, true); });
        rejects_invalid([&] { return stft::inverse_bluestein<Sample>(
                                   transform, invalid); });

        // The Result is user constructible: malformed frequency bins
        // must be rejected before running an inverse FFT.
        auto bad_result = transform;
        bad_result.frames[0][0] = std::complex<T>{bad, T{}};
        rejects_invalid([&] { return stft::inverse_bluestein<Sample>(
                                   bad_result, window); });
        bad_result = transform;
        bad_result.frames[0][0] = std::complex<T>{T{}, bad};
        rejects_invalid([&] { return stft::inverse_bluestein<Sample>(
                                   bad_result, window); });
    }
    // Finite yet nonzero window amplitudes may have unrepresentable
    // squares; inverse normalization must not silently become zero/Inf.
    for (const T bad : {tiny, huge}) {
        auto invalid = window;
        invalid[0] = bad;
        rejects_invalid([&] { return stft::inverse_bluestein<Sample>(
                                   transform, invalid); });
    }
    // Each window-square can be finite while four overlapping
    // contributions overflow the scalar normalization sum.
    {
        const T amplitude = std::sqrt(std::numeric_limits<T>::max() / T{3});
        const std::vector<T> large(4, amplitude);
        stft::Result<Sample> overlap{};
        overlap.frame_size = 4;
        overlap.hop_size = 1;
        overlap.original_size = 7;
        overlap.frames.resize(4, std::vector<std::complex<T>>(4));
        rejects_invalid([&] { return stft::inverse_bluestein<Sample>(
                                   overlap, large); });
    }
    // Genuine zero weights remain legal. A wholly zero window has
    // zero-energy samples and the documented inverse zero convention.
    const std::vector<T> zero_window(window.size(), T{});
    const auto blank = stft::bluestein<Sample>(input, zero_window, 2, true);
    const auto blank_restored = stft::inverse_bluestein<Sample>(
        blank, zero_window);
    assert(blank_restored.size() == input.size());
    for (const auto& sample : blank_restored) assert(std::abs(sample) == T{});
}
}

int main() {
    verify<float>();
    verify<double>();
    verify<std::complex<float>>();
    verify<std::complex<double>>();
}
