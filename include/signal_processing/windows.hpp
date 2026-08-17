#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace signal_processing::windows {

enum class Sampling { symmetric, periodic };

namespace detail {

template <std::floating_point T>
[[nodiscard]] inline T denominator(std::size_t n, Sampling sampling) {
    return static_cast<T>(sampling == Sampling::symmetric ? n - 1 : n);
}

template <std::floating_point T>
[[nodiscard]] inline T normalized_position(std::size_t i, std::size_t n,
                                           Sampling sampling) {
    return T{2} * static_cast<T>(i) / denominator<T>(n, sampling) - T{1};
}

template <std::floating_point T>
[[nodiscard]] inline T phase(std::size_t i, std::size_t n, Sampling sampling) {
    return T{2} * std::numbers::pi_v<T> * static_cast<T>(i) /
           denominator<T>(n, sampling);
}

template <std::floating_point T>
[[nodiscard]] inline T bessel_i0(T x) {
    const T y = x * x / T{4};
    T term{1};
    T sum{1};
    for (std::size_t k = 1; k <= 100; ++k) {
        const T kt = static_cast<T>(k);
        term *= y / (kt * kt);
        sum += term;
        if (std::abs(term) <= std::numeric_limits<T>::epsilon() * std::abs(sum)) break;
    }
    return sum;
}

template <std::floating_point T, typename Function>
[[nodiscard]] inline std::vector<T> generate(std::size_t n, Function&& function) {
    std::vector<T> output(n, T{});
    if (n == 0) return output;
    if (n == 1) {
        output[0] = T{1};
        return output;
    }
    for (std::size_t i = 0; i < n; ++i) output[i] = function(i);
    return output;
}

}  // namespace detail

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> rectangular(std::size_t n) {
    return std::vector<T>(n, T{1});
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> bartlett(std::size_t n,
                                             Sampling sampling = Sampling::symmetric) {
    return detail::generate<T>(n, [=](std::size_t i) {
        return std::max(T{0}, T{1} - std::abs(detail::normalized_position<T>(i, n, sampling)));
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> hann(std::size_t n,
                                         Sampling sampling = Sampling::symmetric) {
    return detail::generate<T>(n, [=](std::size_t i) {
        return T{0.5} - T{0.5} * std::cos(detail::phase<T>(i, n, sampling));
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> hamming(std::size_t n,
                                            Sampling sampling = Sampling::symmetric) {
    constexpr double a0 = 0.54;
    constexpr double a1 = 0.46;
    return detail::generate<T>(n, [=](std::size_t i) {
        return static_cast<T>(a0) - static_cast<T>(a1) *
               std::cos(detail::phase<T>(i, n, sampling));
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> blackman(std::size_t n,
                                             Sampling sampling = Sampling::symmetric) {
    constexpr double a0 = 0.42;
    constexpr double a1 = 0.50;
    constexpr double a2 = 0.08;
    return detail::generate<T>(n, [=](std::size_t i) {
        const T p = detail::phase<T>(i, n, sampling);
        return static_cast<T>(a0) - static_cast<T>(a1) * std::cos(p) +
               static_cast<T>(a2) * std::cos(T{2} * p);
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> blackman_harris(
    std::size_t n, Sampling sampling = Sampling::symmetric) {
    constexpr double a0 = 0.35875;
    constexpr double a1 = 0.48829;
    constexpr double a2 = 0.14128;
    constexpr double a3 = 0.01168;
    return detail::generate<T>(n, [=](std::size_t i) {
        const T p = detail::phase<T>(i, n, sampling);
        return static_cast<T>(a0) - static_cast<T>(a1) * std::cos(p) +
               static_cast<T>(a2) * std::cos(T{2} * p) -
               static_cast<T>(a3) * std::cos(T{3} * p);
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> nuttall(std::size_t n,
                                            Sampling sampling = Sampling::symmetric) {
    constexpr double a0 = 0.355768;
    constexpr double a1 = 0.487396;
    constexpr double a2 = 0.144232;
    constexpr double a3 = 0.012604;
    return detail::generate<T>(n, [=](std::size_t i) {
        const T p = detail::phase<T>(i, n, sampling);
        return static_cast<T>(a0) - static_cast<T>(a1) * std::cos(p) +
               static_cast<T>(a2) * std::cos(T{2} * p) -
               static_cast<T>(a3) * std::cos(T{3} * p);
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> flat_top(std::size_t n,
                                             Sampling sampling = Sampling::symmetric) {
    constexpr double a0 = 0.21557895;
    constexpr double a1 = 0.41663158;
    constexpr double a2 = 0.277263158;
    constexpr double a3 = 0.083578947;
    constexpr double a4 = 0.006947368;
    return detail::generate<T>(n, [=](std::size_t i) {
        const T p = detail::phase<T>(i, n, sampling);
        return static_cast<T>(a0) - static_cast<T>(a1) * std::cos(p) +
               static_cast<T>(a2) * std::cos(T{2} * p) -
               static_cast<T>(a3) * std::cos(T{3} * p) +
               static_cast<T>(a4) * std::cos(T{4} * p);
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> kaiser(std::size_t n, T beta,
                                           Sampling sampling = Sampling::symmetric) {
    if (beta < T{0}) throw std::invalid_argument("Kaiser beta must be nonnegative");
    const T denominator = detail::bessel_i0(beta);
    return detail::generate<T>(n, [=](std::size_t i) {
        const T r = detail::normalized_position<T>(i, n, sampling);
        const T argument = beta * std::sqrt(std::max(T{0}, T{1} - r * r));
        return detail::bessel_i0(argument) / denominator;
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> gaussian(std::size_t n, T sigma,
                                             Sampling sampling = Sampling::symmetric) {
    if (!(sigma > T{0})) throw std::invalid_argument("Gaussian sigma must be positive");
    return detail::generate<T>(n, [=](std::size_t i) {
        const T r = detail::normalized_position<T>(i, n, sampling) / sigma;
        return std::exp(-T{0.5} * r * r);
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> tukey(std::size_t n, T alpha,
                                          Sampling sampling = Sampling::symmetric) {
    if (alpha < T{0} || alpha > T{1}) {
        throw std::invalid_argument("Tukey alpha must lie in [0, 1]");
    }
    if (alpha == T{0}) return rectangular<T>(n);
    if (alpha == T{1}) return hann<T>(n, sampling);

    return detail::generate<T>(n, [=](std::size_t i) {
        const T x = static_cast<T>(i) / detail::denominator<T>(n, sampling);
        if (x < alpha / T{2}) {
            return T{0.5} * (T{1} + std::cos(std::numbers::pi_v<T> *
                                              (T{2} * x / alpha - T{1})));
        }
        if (x <= T{1} - alpha / T{2}) return T{1};
        return T{0.5} * (T{1} + std::cos(std::numbers::pi_v<T> *
                                          (T{2} * x / alpha - T{2} / alpha + T{1})));
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> lanczos(std::size_t n,
                                            Sampling sampling = Sampling::symmetric) {
    return detail::generate<T>(n, [=](std::size_t i) {
        const T r = detail::normalized_position<T>(i, n, sampling);
        if (r == T{0}) return T{1};
        return std::sin(std::numbers::pi_v<T> * r) /
               (std::numbers::pi_v<T> * r);
    });
}

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> welch(std::size_t n,
                                          Sampling sampling = Sampling::symmetric) {
    return detail::generate<T>(n, [=](std::size_t i) {
        const T r = detail::normalized_position<T>(i, n, sampling);
        return std::max(T{0}, T{1} - r * r);
    });
}

}  // namespace signal_processing::windows
