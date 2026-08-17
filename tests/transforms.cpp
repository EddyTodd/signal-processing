#include "signal_processing/signal_processing.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <span>
#include <type_traits>
#include <vector>

namespace {

template <typename T>
[[nodiscard]] T tolerance() {
    if constexpr (std::same_as<T, float>) return 4.0e-3F;
    return 2.0e-9;
}

template <typename T>
void expect_near(std::span<const T> actual, std::span<const T> expected, T multiplier = T{1}) {
    assert(actual.size() == expected.size());
    const T eps = tolerance<T>() * multiplier;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        assert(std::abs(actual[i] - expected[i]) <= eps * (T{1} + std::abs(expected[i])));
    }
}

template <typename T>
[[nodiscard]] std::vector<T> make_signal(std::size_t n) {
    std::vector<T> output(n);
    for (std::size_t i = 0; i < n; ++i) {
        const int value = static_cast<int>((7 * i + 3) % 13) - 6;
        output[i] = static_cast<T>(value) / T{4};
    }
    return output;
}

template <typename T>
void test_dct() {
    namespace dct = signal_processing::dct;
    for (std::size_t n = 1; n <= 9; ++n) {
        const auto input = make_signal<T>(n);
        expect_near<T>(dct::type1_bluestein<T>(input), dct::type1<T>(input), T{6});
        expect_near<T>(dct::type2_bluestein<T>(input), dct::type2<T>(input), T{8});
        expect_near<T>(dct::type3_bluestein<T>(input), dct::type3<T>(input), T{8});
        expect_near<T>(dct::type4_bluestein<T>(input), dct::type4<T>(input), T{12});

        if (n > 1) {
            auto restored1 = dct::type1<T>(dct::type1<T>(input));
            for (auto& value : restored1) value /= T{2} * static_cast<T>(n - 1);
            expect_near<T>(restored1, input, T{3});
        }

        auto restored23 = dct::type3<T>(dct::type2<T>(input));
        for (auto& value : restored23) value /= T{2} * static_cast<T>(n);
        expect_near<T>(restored23, input, T{3});

        auto restored4 = dct::type4<T>(dct::type4<T>(input));
        for (auto& value : restored4) value /= T{2} * static_cast<T>(n);
        expect_near<T>(restored4, input, T{3});
    }
}

template <typename T>
void test_dst() {
    namespace dst = signal_processing::dst;
    for (std::size_t n = 1; n <= 9; ++n) {
        const auto input = make_signal<T>(n);
        expect_near<T>(dst::type1_bluestein<T>(input), dst::type1<T>(input), T{8});
        expect_near<T>(dst::type2_bluestein<T>(input), dst::type2<T>(input), T{8});
        expect_near<T>(dst::type3_bluestein<T>(input), dst::type3<T>(input), T{8});
        expect_near<T>(dst::type4_bluestein<T>(input), dst::type4<T>(input), T{12});

        auto restored1 = dst::type1<T>(dst::type1<T>(input));
        for (auto& value : restored1) value /= T{2} * static_cast<T>(n + 1);
        expect_near<T>(restored1, input, T{3});

        auto restored23 = dst::type3<T>(dst::type2<T>(input));
        for (auto& value : restored23) value /= T{2} * static_cast<T>(n);
        expect_near<T>(restored23, input, T{3});

        auto restored4 = dst::type4<T>(dst::type4<T>(input));
        for (auto& value : restored4) value /= T{2} * static_cast<T>(n);
        expect_near<T>(restored4, input, T{3});
    }
}

template <typename T>
void test_hartley() {
    namespace hartley = signal_processing::hartley;
    for (std::size_t n = 1; n <= 11; ++n) {
        const auto input = make_signal<T>(n);
        const auto direct = hartley::direct<T>(input);
        expect_near<T>(hartley::bluestein<T>(input), direct, T{10});
        auto restored = hartley::direct<T>(direct);
        for (auto& value : restored) value /= static_cast<T>(n);
        expect_near<T>(restored, input, T{4});
    }
}

template <typename T>
void test_walsh_hadamard() {
    namespace wht = signal_processing::walsh_hadamard;
    for (const std::size_t n : {1U, 2U, 4U, 8U, 16U}) {
        const auto input = make_signal<T>(n);
        const auto direct = wht::direct<T>(input);
        expect_near<T>(wht::fast<T>(input), direct);
        auto restored = wht::fast<T>(direct);
        for (auto& value : restored) value /= static_cast<T>(n);
        expect_near<T>(restored, input);
    }
}

template <typename T>
void test_mdct() {
    namespace mdct = signal_processing::mdct;
    constexpr std::size_t n = 8;
    constexpr std::size_t selected = 3;
    std::vector<T> basis(2 * n);
    for (std::size_t j = 0; j < 2 * n; ++j) {
        const T angle = std::numbers::pi_v<T> / static_cast<T>(n) *
                        (static_cast<T>(j) + T{0.5} + static_cast<T>(n) / T{2}) *
                        (static_cast<T>(selected) + T{0.5});
        basis[j] = std::cos(angle);
    }

    const auto spectrum = mdct::direct<T>(basis);
    for (std::size_t k = 0; k < n; ++k) {
        const T expected = k == selected ? static_cast<T>(n) : T{0};
        assert(std::abs(spectrum[k] - expected) <= tolerance<T>() * T{8} *
                                                      (T{1} + std::abs(expected)));
    }

    std::vector<T> impulse(n, T{});
    impulse[selected] = T{1};
    const auto inverse = mdct::inverse_direct<T>(impulse);
    for (std::size_t j = 0; j < 2 * n; ++j) {
        const T expected = T{2} / static_cast<T>(n) * basis[j];
        assert(std::abs(inverse[j] - expected) <= tolerance<T>() * T{4});
    }
}

template <typename T>
void test_windows() {
    namespace windows = signal_processing::windows;
    const auto hann = windows::hann<T>(5);
    assert(std::abs(hann.front()) <= tolerance<T>());
    assert(std::abs(hann[2] - T{1}) <= tolerance<T>());
    assert(std::abs(hann.back()) <= tolerance<T>());

    const auto periodic_hann = windows::hann<T>(5, windows::Sampling::periodic);
    assert(std::abs(periodic_hann.back()) > tolerance<T>());

    const auto bartlett = windows::bartlett<T>(5);
    assert(std::abs(bartlett.front()) <= tolerance<T>());
    assert(std::abs(bartlett[2] - T{1}) <= tolerance<T>());
    assert(std::abs(bartlett.back()) <= tolerance<T>());

    const auto rectangular = windows::rectangular<T>(7);
    const auto kaiser_zero = windows::kaiser<T>(7, T{0});
    const auto tukey_zero = windows::tukey<T>(7, T{0});
    expect_near<T>(kaiser_zero, rectangular);
    expect_near<T>(tukey_zero, rectangular);
    expect_near<T>(windows::tukey<T>(7, T{1}), windows::hann<T>(7));

    const auto gaussian = windows::gaussian<T>(7, T{0.4});
    const auto kaiser = windows::kaiser<T>(7, T{8});
    for (std::size_t i = 0; i < 7; ++i) {
        assert(std::abs(gaussian[i] - gaussian[6 - i]) <= tolerance<T>());
        assert(std::abs(kaiser[i] - kaiser[6 - i]) <= tolerance<T>());
    }

    const auto flat_top = windows::flat_top<T>(5);
    assert(std::abs(flat_top[2] - T{1}) <= tolerance<T>() * T{2});
    const auto lanczos = windows::lanczos<T>(5);
    assert(std::abs(lanczos.front()) <= tolerance<T>());
    assert(std::abs(lanczos[2] - T{1}) <= tolerance<T>());
    const auto welch = windows::welch<T>(5);
    assert(std::abs(welch.front()) <= tolerance<T>());
    assert(std::abs(welch[2] - T{1}) <= tolerance<T>());

    assert(windows::blackman_harris<T>(1)[0] == T{1});
    assert(windows::nuttall<T>(1)[0] == T{1});
}

}  // namespace

int main() {
    test_dct<float>();
    test_dct<double>();
    test_dst<float>();
    test_dst<double>();
    test_hartley<float>();
    test_hartley<double>();
    test_walsh_hadamard<float>();
    test_walsh_hadamard<double>();
    test_mdct<float>();
    test_mdct<double>();
    test_windows<float>();
    test_windows<double>();
}
