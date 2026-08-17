#pragma once

#include "signal_processing/iir.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace signal_processing::iir::design {

template <Scalar T>
using Complex = std::complex<T>;

template <Scalar T>
struct Zpk {
    std::vector<Complex<T>> zeros;
    std::vector<Complex<T>> poles;
    T gain{1};
};

enum class Response { lowpass, highpass, bandpass, bandstop };
enum class BesselNormalization { phase, delay, magnitude_3db };

template <Scalar T>
struct SosDesign {
    std::vector<BiquadCoefficients<T>> sections;
    T gain{1};
};

namespace detail {

template <Scalar T>
[[nodiscard]] T pow10m1(T x) {
    return std::expm1(std::numbers::ln10_v<T> * x);
}

template <Scalar T>
[[nodiscard]] Complex<T> product_negated(std::span<const Complex<T>> roots) {
    Complex<T> result{1, 0};
    for (const auto& root : roots) result *= -root;
    return result;
}

template <Scalar T>
[[nodiscard]] std::size_t relative_degree(const Zpk<T>& zpk) {
    if (zpk.zeros.size() > zpk.poles.size())
        throw std::invalid_argument("ZPK must have at least as many poles as zeros");
    return zpk.poles.size() - zpk.zeros.size();
}

[[nodiscard]] inline std::size_t twice_count(std::size_t count, const char* message) {
    if (count > std::numeric_limits<std::size_t>::max() / 2)
        throw std::length_error(message);
    return 2 * count;
}

[[nodiscard]] inline std::size_t twice_plus(std::size_t count, std::size_t extra,
                                             const char* message) {
    if (extra > std::numeric_limits<std::size_t>::max() ||
        count > (std::numeric_limits<std::size_t>::max() - extra) / 2)
        throw std::length_error(message);
    return 2 * count + extra;
}

[[nodiscard]] inline std::size_t twice_sum(std::size_t first, std::size_t second,
                                            const char* message) {
    if (first > std::numeric_limits<std::size_t>::max() - second)
        throw std::length_error(message);
    return twice_count(first + second, message);
}

template <Scalar T>
[[nodiscard]] T carlson_rf(T x, T y, T z) {
    if (x < T{} || y < T{} || z < T{} || (x + y == T{}) ||
        (x + z == T{}) || (y + z == T{}))
        throw std::invalid_argument("Carlson RF requires nonnegative arguments with at most one zero");
    constexpr T c1 = static_cast<T>(1.0 / 24.0);
    constexpr T c2 = static_cast<T>(0.1);
    constexpr T c3 = static_cast<T>(3.0 / 44.0);
    constexpr T c4 = static_cast<T>(1.0 / 14.0);
    const T tolerance = std::pow(T{4} * std::numeric_limits<T>::epsilon(), T{1} / T{6});
    T average{}, dx{}, dy{}, dz{};
    for (int iteration = 0; iteration < 64; ++iteration) {
        const T sx = std::sqrt(x);
        const T sy = std::sqrt(y);
        const T sz = std::sqrt(z);
        const T lambda = sx * (sy + sz) + sy * sz;
        x = (x + lambda) / T{4};
        y = (y + lambda) / T{4};
        z = (z + lambda) / T{4};
        average = (x + y + z) / T{3};
        dx = (average - x) / average;
        dy = (average - y) / average;
        dz = (average - z) / average;
        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) < tolerance) {
            const T e2 = dx * dy - dz * dz;
            const T e3 = dx * dy * dz;
            return (T{1} + (c1 * e2 - c2 - c3 * e3) * e2 + c4 * e3) /
                   std::sqrt(average);
        }
    }
    throw std::runtime_error("Carlson RF did not converge");
}

template <Scalar T>
[[nodiscard]] T incomplete_elliptic_f(T phi, T m) {
    const T s = std::sin(phi);
    const T c = std::cos(phi);
    return s * carlson_rf<T>(c * c, T{1} - m * s * s, T{1});
}

template <Scalar T>
[[nodiscard]] T complete_elliptic_k(T m) {
    if (!(m >= T{} && m < T{1}))
        throw std::invalid_argument("elliptic parameter must be in [0, 1)");
    return carlson_rf<T>(T{}, T{1} - m, T{1});
}

template <Scalar T>
struct JacobiValues {
    T sn{};
    T cn{};
    T dn{};
};

template <Scalar T>
[[nodiscard]] JacobiValues<T> jacobi_real(T u, T m) {
    if (m < T{} || m >= T{1})
        throw std::invalid_argument("Jacobi parameter must be in [0, 1)");
    if (m == T{}) return {std::sin(u), std::cos(u), T{1}};
    const T capk = complete_elliptic_k<T>(m);
    if (std::abs(u) > capk)
        throw std::invalid_argument("Jacobi helper supports |u| <= K(m)");
    const T sign = u < T{} ? T{-1} : T{1};
    const T target = std::abs(u);
    T lo{};
    T hi = std::numbers::pi_v<T> / T{2};
    for (int iteration = 0; iteration < 96; ++iteration) {
        const T mid = (lo + hi) / T{2};
        if (incomplete_elliptic_f<T>(mid, m) < target) lo = mid;
        else hi = mid;
    }
    const T phi = sign * (lo + hi) / T{2};
    const T sn = std::sin(phi);
    const T cn = std::cos(phi);
    return {sn, cn, std::sqrt(std::max(T{}, T{1} - m * sn * sn))};
}

template <Scalar T>
[[nodiscard]] T inverse_sc_complement(T w, T m1) {
    if (!(w >= T{}) || !(m1 > T{} && m1 < T{1}))
        throw std::invalid_argument("inverse sc requires w >= 0 and m in (0,1)");
    const T m = T{1} - m1;
    T lo{};
    T hi = complete_elliptic_k<T>(m) *
           (T{1} - T{64} * std::numeric_limits<T>::epsilon());
    for (int iteration = 0; iteration < 96; ++iteration) {
        const T mid = (lo + hi) / T{2};
        const auto values = jacobi_real<T>(mid, m);
        const T sc = values.sn / values.cn;
        if (sc < w) lo = mid;
        else hi = mid;
    }
    return (lo + hi) / T{2};
}

template <Scalar T>
[[nodiscard]] T elliptic_degree(std::size_t n, T m1) {
    const T k1 = complete_elliptic_k<T>(m1);
    const T k1p = complete_elliptic_k<T>(T{1} - m1);
    const T q1 = std::exp(-std::numbers::pi_v<T> * k1p / k1);
    const T q = std::pow(q1, T{1} / static_cast<T>(n));
    T numerator{};
    for (int j = 0; j <= 7; ++j) numerator += std::pow(q, static_cast<T>(j * (j + 1)));
    T denominator{1};
    for (int j = 1; j <= 8; ++j) denominator += T{2} * std::pow(q, static_cast<T>(j * j));
    const T ratio = numerator / denominator;
    return T{16} * q * ratio * ratio * ratio * ratio;
}

template <Scalar T>
[[nodiscard]] Complex<T> polynomial_value(std::span<const T> coefficients, Complex<T> z) {
    Complex<T> result{};
    for (std::size_t i = coefficients.size(); i-- > 0;) result = result * z + coefficients[i];
    return result;
}

template <Scalar T>
[[nodiscard]] Complex<T> polynomial_derivative(std::span<const T> coefficients,
                                                Complex<T> z) {
    Complex<T> result{};
    for (std::size_t i = coefficients.size(); i-- > 1;)
        result = result * z + static_cast<T>(i) * coefficients[i];
    return result;
}

template <Scalar T>
[[nodiscard]] std::vector<Complex<T>> scaled_polynomial_roots_aberth(
    std::span<const T> coefficients) {
    if (coefficients.empty())
        throw std::invalid_argument("scaled root finder requires a nonempty polynomial");
    const std::size_t n = coefficients.size() - 1;
    if (n == 0 || coefficients.front() == T{} || coefficients.back() == T{})
        throw std::invalid_argument("scaled root finder requires nonzero constant and leading coefficients");
    const T scale = std::pow(std::abs(coefficients.front() / coefficients.back()),
                             T{1} / static_cast<T>(n));
    std::vector<T> scaled(coefficients.size());
    T power{1};
    for (std::size_t k = 0; k < coefficients.size(); ++k) {
        scaled[k] = coefficients[k] * power / coefficients.front();
        power *= scale;
    }
    std::vector<Complex<T>> roots(n);
    for (std::size_t i = 0; i < n; ++i) {
        const T angle = T{2} * std::numbers::pi_v<T> *
                        (static_cast<T>(i) + T{0.5}) / static_cast<T>(n);
        roots[i] = std::polar(T{1}, angle);
    }
    const T tolerance = T{4096} * std::numeric_limits<T>::epsilon();
    for (int iteration = 0; iteration < 256; ++iteration) {
        for (std::size_t i = 0; i < n; ++i) {
            const auto f = polynomial_value<T>(scaled, roots[i]);
            const auto fp = polynomial_derivative<T>(scaled, roots[i]);
            Complex<T> sum{};
            for (std::size_t j = 0; j < n; ++j)
                if (j != i) sum += T{1} / (roots[i] - roots[j]);
            const auto denominator = fp - f * sum;
            if (std::abs(denominator) != T{}) roots[i] -= f / denominator;
        }
        T maximum_residual{};
        for (const auto& root : roots) {
            const T numerator = std::abs(polynomial_value<T>(scaled, root));
            T denominator{};
            T magnitude_power{1};
            const T magnitude = std::abs(root);
            for (const T coefficient : scaled) {
                denominator += std::abs(coefficient) * magnitude_power;
                magnitude_power *= magnitude;
            }
            maximum_residual = std::max(maximum_residual,
                                        numerator / std::max(denominator, T{1}));
        }
        if (maximum_residual <= tolerance) {
            for (auto& root : roots) root *= scale;
            return roots;
        }
    }
    throw std::runtime_error("Bessel polynomial root iteration did not converge");
}

template <Scalar T>
[[nodiscard]] std::vector<T> reverse_bessel_coefficients(std::size_t n) {
    if (n == std::numeric_limits<std::size_t>::max())
        throw std::length_error("Bessel coefficient count overflows size_t");
    if (n == 0) return {T{1}};
    std::vector<T> coefficients(n + 1);
    T a0{1};
    for (std::size_t j = 1; j <= n; ++j)
        a0 *= (static_cast<T>(n) + static_cast<T>(j)) / T{2};
    coefficients[0] = a0;
    for (std::size_t k = 0; k < n; ++k) {
        coefficients[k + 1] = coefficients[k] * T{2} * static_cast<T>(n - k) /
                              (static_cast<T>(k + 1) *
                               (T{2} * static_cast<T>(n) - static_cast<T>(k)));
    }
    return coefficients;
}

template <Scalar T>
[[nodiscard]] Complex<T> analog_response(const Zpk<T>& zpk, Complex<T> s) {
    Complex<T> numerator{zpk.gain, 0};
    Complex<T> denominator{1, 0};
    for (const auto& zero : zpk.zeros) numerator *= s - zero;
    for (const auto& pole : zpk.poles) denominator *= s - pole;
    return numerator / denominator;
}

template <Scalar T>
[[nodiscard]] std::vector<Complex<T>> polynomial_from_roots(
    std::span<const Complex<T>> roots) {
    std::vector<Complex<T>> coefficients{Complex<T>{1, 0}};
    for (const auto& root : roots) {
        coefficients.push_back({});
        for (std::size_t i = coefficients.size() - 1; i > 0; --i)
            coefficients[i] = coefficients[i] - root * coefficients[i - 1];
    }
    return coefficients;
}

template <Scalar T>
[[nodiscard]] std::vector<std::array<Complex<T>, 2>> pair_roots(
    std::vector<Complex<T>> roots) {
    const T tolerance = T{256} * std::numeric_limits<T>::epsilon();
    std::vector<std::array<Complex<T>, 2>> pairs;

    while (!roots.empty()) {
        const auto first = roots.back();
        roots.pop_back();
        const bool first_is_real =
            std::abs(first.imag()) <= tolerance * (T{1} + std::abs(first.real()));

        if (first_is_real) {
            const auto it = std::find_if(roots.begin(), roots.end(), [&](const auto& root) {
                return std::abs(root.imag()) <= tolerance * (T{1} + std::abs(root.real()));
            });
            if (it == roots.end()) {
                pairs.push_back({Complex<T>{first.real(), T{}}, Complex<T>{T{}, T{}}});
            } else {
                const auto second = *it;
                roots.erase(it);
                pairs.push_back({Complex<T>{first.real(), T{}}, Complex<T>{second.real(), T{}}});
            }
            continue;
        }

        const auto target = std::conj(first);
        const auto it = std::min_element(roots.begin(), roots.end(), [&](const auto& lhs,
                                                                         const auto& rhs) {
            return std::abs(lhs - target) < std::abs(rhs - target);
        });
        if (it == roots.end() ||
            std::abs(*it - target) > static_cast<T>(1e-3) * (T{1} + std::abs(target)))
            throw std::runtime_error("could not pair conjugate filter roots");
        const auto second = *it;
        roots.erase(it);
        pairs.push_back({first, second});
    }
    return pairs;
}

template <Scalar T>
[[nodiscard]] BiquadCoefficients<T> section_from_pairs(
    const std::array<Complex<T>, 2>& zeros, const std::array<Complex<T>, 2>& poles) {
    const auto zero_sum = zeros[0] + zeros[1];
    const auto zero_product = zeros[0] * zeros[1];
    const auto pole_sum = poles[0] + poles[1];
    const auto pole_product = poles[0] * poles[1];
    return {T{1}, -zero_sum.real(), zero_product.real(), T{1}, -pole_sum.real(),
            pole_product.real()};
}

}  // namespace detail

// Butterworth filter.
// https://en.wikipedia.org/wiki/Butterworth_filter
template <Scalar T>
[[nodiscard]] Zpk<T> butterworth_prototype(std::size_t order) {
    if (order == 0) return {};
    Zpk<T> result;
    result.poles.reserve(order);
    for (std::size_t k = 0; k < order; ++k) {
        const T m = T{1} + T{2} * static_cast<T>(k) - static_cast<T>(order);
        const T angle = std::numbers::pi_v<T> * m / (T{2} * static_cast<T>(order));
        result.poles.push_back(-std::exp(Complex<T>{T{}, angle}));
    }
    result.gain = T{1};
    return result;
}

// Chebyshev filter.
// https://en.wikipedia.org/wiki/Chebyshev_filter
template <Scalar T>
[[nodiscard]] Zpk<T> chebyshev1_prototype(std::size_t order, T ripple_db) {
    if (!(ripple_db > T{}) || !std::isfinite(ripple_db))
        throw std::invalid_argument("Chebyshev-I ripple must be finite and positive");
    if (order == 0) return {{}, {}, std::pow(T{10}, -ripple_db / T{20})};

    const T epsilon = std::sqrt(detail::pow10m1<T>(ripple_db / T{10}));
    const T mu = std::asinh(T{1} / epsilon) / static_cast<T>(order);
    Zpk<T> result;
    result.poles.reserve(order);
    for (std::size_t k = 0; k < order; ++k) {
        const T m = T{1} + T{2} * static_cast<T>(k) - static_cast<T>(order);
        const T theta = std::numbers::pi_v<T> * m / (T{2} * static_cast<T>(order));
        result.poles.push_back(-std::sinh(Complex<T>{mu, theta}));
    }
    result.gain = detail::product_negated<T>(result.poles).real();
    if ((order & 1U) == 0U) result.gain /= std::sqrt(T{1} + epsilon * epsilon);
    return result;
}

template <Scalar T>
[[nodiscard]] Zpk<T> chebyshev2_prototype(std::size_t order, T stop_db) {
    if (!(stop_db > T{}) || !std::isfinite(stop_db))
        throw std::invalid_argument("Chebyshev-II attenuation must be finite and positive");
    if (order == 0) return {};

    const T de = T{1} / std::sqrt(detail::pow10m1<T>(stop_db / T{10}));
    const T mu = std::asinh(T{1} / de) / static_cast<T>(order);
    Zpk<T> result;
    result.poles.reserve(order);
    result.zeros.reserve(order);
    for (std::size_t k = 0; k < order; ++k) {
        const T m = T{1} + T{2} * static_cast<T>(k) - static_cast<T>(order);
        const T theta = std::numbers::pi_v<T> * m / (T{2} * static_cast<T>(order));
        if (!((order & 1U) != 0U && m == T{}))
            result.zeros.emplace_back(T{}, T{1} / std::sin(theta));
        result.poles.push_back(-T{1} / std::sinh(Complex<T>{mu, theta}));
    }
    result.gain = (detail::product_negated<T>(result.poles) /
                   detail::product_negated<T>(result.zeros)).real();
    return result;
}

// Elliptic filter.
// https://en.wikipedia.org/wiki/Elliptic_filter
template <Scalar T>
[[nodiscard]] Zpk<T> elliptic_prototype(std::size_t order, T ripple_db, T stop_db) {
    if (!(ripple_db > T{} && stop_db > ripple_db) ||
        !std::isfinite(ripple_db) || !std::isfinite(stop_db))
        throw std::invalid_argument("elliptic design requires finite 0 < ripple_db < stop_db");
    if (order == 0) return {{}, {}, std::pow(T{10}, -ripple_db / T{20})};
    if (order == 1) {
        const T pole = -std::sqrt(T{1} / detail::pow10m1<T>(ripple_db / T{10}));
        return {{}, {Complex<T>{pole, T{}}}, -pole};
    }

    const T epsilon_squared = detail::pow10m1<T>(ripple_db / T{10});
    const T complementary_modulus_squared =
        epsilon_squared / detail::pow10m1<T>(stop_db / T{10});
    if (!(complementary_modulus_squared > T{} && complementary_modulus_squared < T{1}))
        throw std::invalid_argument("elliptic specifications are not realizable");

    const T epsilon = std::sqrt(epsilon_squared);
    const T modulus_squared = detail::elliptic_degree<T>(order, complementary_modulus_squared);
    const T complete_k = detail::complete_elliptic_k<T>(modulus_squared);

    std::vector<T> indices;
    for (std::size_t j = (order & 1U) != 0U ? 0U : 1U; j < order; j += 2)
        indices.push_back(static_cast<T>(j));

    struct JacobiTriple { T sn; T cn; T dn; };
    std::vector<JacobiTriple> values;
    values.reserve(indices.size());
    for (const T j : indices) {
        const auto value =
            detail::jacobi_real<T>(j * complete_k / static_cast<T>(order), modulus_squared);
        values.push_back({value.sn, value.cn, value.dn});
    }

    Zpk<T> result;
    for (const auto& value : values) {
        if (std::abs(value.sn) > T{32} * std::numeric_limits<T>::epsilon()) {
            const T imaginary = T{1} / (std::sqrt(modulus_squared) * value.sn);
            result.zeros.emplace_back(T{}, imaginary);
            result.zeros.emplace_back(T{}, -imaginary);
        }
    }

    const T inverse_sc =
        detail::inverse_sc_complement<T>(T{1} / epsilon, complementary_modulus_squared);
    const T v0 = complete_k * inverse_sc /
                 (static_cast<T>(order) *
                  detail::complete_elliptic_k<T>(complementary_modulus_squared));
    const auto v = detail::jacobi_real<T>(v0, T{1} - modulus_squared);

    std::vector<Complex<T>> base_poles;
    base_poles.reserve(values.size());
    for (const auto& value : values) {
        const T denominator = T{1} - (value.dn * v.sn) * (value.dn * v.sn);
        base_poles.emplace_back(-(value.cn * value.dn * v.sn * v.cn) / denominator,
                                -(value.sn * v.dn) / denominator);
    }

    if ((order & 1U) != 0U) {
        result.poles.emplace_back(base_poles.front().real(), T{});
        for (std::size_t i = 1; i < base_poles.size(); ++i) {
            result.poles.push_back(base_poles[i]);
            result.poles.push_back(std::conj(base_poles[i]));
        }
    } else {
        for (const auto& pole : base_poles) {
            result.poles.push_back(pole);
            result.poles.push_back(std::conj(pole));
        }
    }

    result.gain = (detail::product_negated<T>(result.poles) /
                   detail::product_negated<T>(result.zeros)).real();
    if ((order & 1U) == 0U) result.gain /= std::sqrt(T{1} + epsilon_squared);
    return result;
}

// Bessel filter.
// https://en.wikipedia.org/wiki/Bessel_filter
template <Scalar T>
[[nodiscard]] Zpk<T> bessel_prototype(
    std::size_t order, BesselNormalization normalization = BesselNormalization::phase) {
    if (order == 0) return {};

    const auto coefficients = detail::reverse_bessel_coefficients<T>(order);
    auto poles = detail::scaled_polynomial_roots_aberth<T>(coefficients);
    std::sort(poles.begin(), poles.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.imag() < rhs.imag();
    });

    for (std::size_t i = 0; i < poles.size() / 2; ++i) {
        const std::size_t j = poles.size() - 1 - i;
        const T real = (poles[i].real() + poles[j].real()) / T{2};
        const T imaginary =
            (std::abs(poles[i].imag()) + std::abs(poles[j].imag())) / T{2};
        poles[i] = {real, -imaginary};
        poles[j] = {real, imaginary};
    }
    if ((poles.size() & 1U) != 0U)
        poles[poles.size() / 2] = {poles[poles.size() / 2].real(), T{}};
    if (std::any_of(poles.begin(), poles.end(),
                    [](const auto& pole) { return pole.real() >= T{}; }))
        throw std::runtime_error("Bessel root solver produced a non-left-half-plane pole");

    const T a0 = coefficients.front();
    T gain = a0;
    switch (normalization) {
        case BesselNormalization::phase: {
            const T scale = std::pow(a0, -T{1} / static_cast<T>(order));
            for (auto& pole : poles) pole *= scale;
            gain = T{1};
            break;
        }
        case BesselNormalization::delay:
            break;
        case BesselNormalization::magnitude_3db: {
            const Zpk<T> natural{{}, poles, a0};
            const T target = T{1} / std::sqrt(T{2});
            T low{};
            T high{1};
            while (std::abs(detail::analog_response<T>(natural, {T{}, high})) > target) {
                if (high > std::numeric_limits<T>::max() / T{2})
                    throw std::runtime_error("Bessel -3 dB normalization failed to bracket the edge");
                high *= T{2};
            }
            for (int iteration = 0; iteration < 96; ++iteration) {
                const T middle = (low + high) / T{2};
                if (std::abs(detail::analog_response<T>(natural, {T{}, middle})) > target)
                    low = middle;
                else
                    high = middle;
            }
            const T w3 = (low + high) / T{2};
            for (auto& pole : poles) pole /= w3;
            gain = a0 / std::pow(w3, static_cast<T>(order));
            break;
        }
        default:
            throw std::invalid_argument("invalid Bessel normalization");
    }
    return {{}, std::move(poles), gain};
}

template <Scalar T>
[[nodiscard]] Zpk<T> lowpass_transform(Zpk<T> zpk, T omega) {
    if (!(omega > T{}) || !std::isfinite(omega))
        throw std::invalid_argument("lowpass frequency must be finite and positive");
    const auto degree = detail::relative_degree<T>(zpk);
    for (auto& zero : zpk.zeros) zero *= omega;
    for (auto& pole : zpk.poles) pole *= omega;
    zpk.gain *= std::pow(omega, static_cast<T>(degree));
    return zpk;
}

template <Scalar T>
[[nodiscard]] Zpk<T> highpass_transform(Zpk<T> zpk, T omega) {
    if (!(omega > T{}) || !std::isfinite(omega))
        throw std::invalid_argument("highpass frequency must be finite and positive");
    const auto degree = detail::relative_degree<T>(zpk);
    const auto old_zeros = zpk.zeros;
    const auto old_poles = zpk.poles;
    for (auto& zero : zpk.zeros) zero = omega / zero;
    for (auto& pole : zpk.poles) pole = omega / pole;
    zpk.zeros.insert(zpk.zeros.end(), degree, Complex<T>{T{}, T{}});
    zpk.gain *=
        (detail::product_negated<T>(old_zeros) / detail::product_negated<T>(old_poles)).real();
    return zpk;
}

template <Scalar T>
[[nodiscard]] Zpk<T> bandpass_transform(Zpk<T> zpk, T omega0, T bandwidth) {
    if (!(omega0 > T{} && bandwidth > T{}) ||
        !std::isfinite(omega0) || !std::isfinite(bandwidth))
        throw std::invalid_argument("bandpass frequencies must be finite and positive");
    const auto degree = detail::relative_degree<T>(zpk);
    std::vector<Complex<T>> zeros;
    std::vector<Complex<T>> poles;
    zeros.reserve(detail::twice_plus(zpk.zeros.size(), degree,
                                      "bandpass zero count overflows size_t"));
    poles.reserve(detail::twice_count(zpk.poles.size(),
                                      "bandpass pole count overflows size_t"));

    for (const auto root : zpk.zeros) {
        const auto q = root * (bandwidth / T{2});
        const auto d = std::sqrt(q * q - Complex<T>{omega0 * omega0, T{}});
        zeros.push_back(q + d);
        zeros.push_back(q - d);
    }
    for (const auto root : zpk.poles) {
        const auto q = root * (bandwidth / T{2});
        const auto d = std::sqrt(q * q - Complex<T>{omega0 * omega0, T{}});
        poles.push_back(q + d);
        poles.push_back(q - d);
    }
    zeros.insert(zeros.end(), degree, Complex<T>{T{}, T{}});
    zpk.zeros = std::move(zeros);
    zpk.poles = std::move(poles);
    zpk.gain *= std::pow(bandwidth, static_cast<T>(degree));
    return zpk;
}

template <Scalar T>
[[nodiscard]] Zpk<T> bandstop_transform(Zpk<T> zpk, T omega0, T bandwidth) {
    if (!(omega0 > T{} && bandwidth > T{}) ||
        !std::isfinite(omega0) || !std::isfinite(bandwidth))
        throw std::invalid_argument("bandstop frequencies must be finite and positive");
    const auto degree = detail::relative_degree<T>(zpk);
    const auto old_zeros = zpk.zeros;
    const auto old_poles = zpk.poles;
    std::vector<Complex<T>> zeros;
    std::vector<Complex<T>> poles;
    zeros.reserve(detail::twice_sum(zpk.zeros.size(), degree,
                                     "bandstop zero count overflows size_t"));
    poles.reserve(detail::twice_count(zpk.poles.size(),
                                      "bandstop pole count overflows size_t"));

    for (const auto root : zpk.zeros) {
        const auto q = (bandwidth / T{2}) / root;
        const auto d = std::sqrt(q * q - Complex<T>{omega0 * omega0, T{}});
        zeros.push_back(q + d);
        zeros.push_back(q - d);
    }
    for (const auto root : zpk.poles) {
        const auto q = (bandwidth / T{2}) / root;
        const auto d = std::sqrt(q * q - Complex<T>{omega0 * omega0, T{}});
        poles.push_back(q + d);
        poles.push_back(q - d);
    }
    for (std::size_t i = 0; i < degree; ++i) {
        zeros.emplace_back(T{}, omega0);
        zeros.emplace_back(T{}, -omega0);
    }
    zpk.zeros = std::move(zeros);
    zpk.poles = std::move(poles);
    zpk.gain *=
        (detail::product_negated<T>(old_zeros) / detail::product_negated<T>(old_poles)).real();
    return zpk;
}

// Bilinear transform.
// https://en.wikipedia.org/wiki/Bilinear_transform
template <Scalar T>
[[nodiscard]] Zpk<T> bilinear_transform(Zpk<T> zpk, T sample_rate = T{1}) {
    if (!(sample_rate > T{}) || !std::isfinite(sample_rate))
        throw std::invalid_argument("sample rate must be finite and positive");
    const auto degree = detail::relative_degree<T>(zpk);
    const T twice_sample_rate = T{2} * sample_rate;
    if (!std::isfinite(twice_sample_rate))
        throw std::overflow_error("twice the sample rate is not representable");
    const auto old_zeros = zpk.zeros;
    const auto old_poles = zpk.poles;

    for (auto& zero : zpk.zeros) zero = (twice_sample_rate + zero) / (twice_sample_rate - zero);
    for (auto& pole : zpk.poles) pole = (twice_sample_rate + pole) / (twice_sample_rate - pole);
    zpk.zeros.insert(zpk.zeros.end(), degree, Complex<T>{T{-1}, T{}});

    Complex<T> ratio{1, 0};
    for (const auto zero : old_zeros) ratio *= twice_sample_rate - zero;
    for (const auto pole : old_poles) ratio /= twice_sample_rate - pole;
    zpk.gain *= ratio.real();
    return zpk;
}

template <Scalar T>
[[nodiscard]] Zpk<T> digital_from_prototype(Zpk<T> prototype, Response response, T first,
                                             T second = T{}) {
    if (!(first > T{} && first < T{1}))
        throw std::invalid_argument("digital critical frequencies are normalized to Nyquist and must be in (0,1)");
    const auto prewarp = [](T frequency) {
        return T{2} * std::tan(std::numbers::pi_v<T> * frequency / T{2});
    };

    switch (response) {
        case Response::lowpass:
            prototype = lowpass_transform<T>(std::move(prototype), prewarp(first));
            break;
        case Response::highpass:
            prototype = highpass_transform<T>(std::move(prototype), prewarp(first));
            break;
        case Response::bandpass:
        case Response::bandstop: {
            if (!(second > first && second < T{1}))
                throw std::invalid_argument("band edges must satisfy 0 < first < second < 1");
            const T w1 = prewarp(first);
            const T w2 = prewarp(second);
            const T center = std::sqrt(w1 * w2);
            const T bandwidth = w2 - w1;
            if (response == Response::bandpass)
                prototype = bandpass_transform<T>(std::move(prototype), center, bandwidth);
            else
                prototype = bandstop_transform<T>(std::move(prototype), center, bandwidth);
            break;
        }
        default:
            throw std::invalid_argument("invalid IIR response type");
    }
    return bilinear_transform<T>(std::move(prototype));
}

template <Scalar T>
[[nodiscard]] Complex<T> analog_frequency_response(const Zpk<T>& zpk, T angular_frequency) {
    if (!(angular_frequency >= T{}) || !std::isfinite(angular_frequency))
        throw std::invalid_argument("analog frequency must be finite and nonnegative");
    return detail::analog_response<T>(zpk, Complex<T>{T{}, angular_frequency});
}

template <Scalar T>
[[nodiscard]] bool analog_stable(const Zpk<T>& zpk) noexcept {
    return std::all_of(zpk.poles.begin(), zpk.poles.end(),
                       [](const auto& pole) { return pole.real() < T{}; });
}

template <Scalar T>
[[nodiscard]] Zpk<T> butterworth(std::size_t order, Response response, T first,
                                 T second = T{}) {
    return digital_from_prototype<T>(butterworth_prototype<T>(order), response, first, second);
}

template <Scalar T>
[[nodiscard]] Zpk<T> chebyshev1(std::size_t order, T ripple_db, Response response, T first,
                                T second = T{}) {
    return digital_from_prototype<T>(chebyshev1_prototype<T>(order, ripple_db), response, first,
                                     second);
}

template <Scalar T>
[[nodiscard]] Zpk<T> chebyshev2(std::size_t order, T stop_db, Response response, T first,
                                T second = T{}) {
    return digital_from_prototype<T>(chebyshev2_prototype<T>(order, stop_db), response, first,
                                     second);
}

template <Scalar T>
[[nodiscard]] Zpk<T> elliptic(std::size_t order, T ripple_db, T stop_db, Response response,
                              T first, T second = T{}) {
    return digital_from_prototype<T>(elliptic_prototype<T>(order, ripple_db, stop_db), response,
                                     first, second);
}

template <Scalar T>
[[nodiscard]] Zpk<T> bessel(
    std::size_t order, Response response, T first, T second = T{},
    BesselNormalization normalization = BesselNormalization::phase) {
    return digital_from_prototype<T>(bessel_prototype<T>(order, normalization), response, first,
                                     second);
}

template <Scalar T>
[[nodiscard]] Complex<T> frequency_response(const Zpk<T>& zpk, T normalized_frequency) {
    if (!(normalized_frequency >= T{} && normalized_frequency <= T{1}))
        throw std::invalid_argument("frequency must be normalized to Nyquist in [0,1]");
    const Complex<T> z =
        std::exp(Complex<T>{T{}, std::numbers::pi_v<T> * normalized_frequency});
    Complex<T> numerator{zpk.gain, T{}};
    Complex<T> denominator{T{1}, T{}};
    for (const auto root : zpk.zeros) numerator *= z - root;
    for (const auto root : zpk.poles) denominator *= z - root;
    return numerator / denominator;
}

template <Scalar T>
[[nodiscard]] bool stable(const Zpk<T>& zpk) noexcept {
    return std::all_of(zpk.poles.begin(), zpk.poles.end(),
                       [](const auto pole) { return std::abs(pole) < T{1}; });
}

template <Scalar T>
[[nodiscard]] Coefficients<T> transfer_function(const Zpk<T>& zpk) {
    if (zpk.zeros.size() != zpk.poles.size())
        throw std::invalid_argument("digital transfer conversion requires equal zero/pole degree");
    const auto numerator_complex = detail::polynomial_from_roots<T>(zpk.zeros);
    const auto denominator_complex = detail::polynomial_from_roots<T>(zpk.poles);
    Coefficients<T> output;
    output.numerator.resize(numerator_complex.size());
    output.denominator.resize(denominator_complex.size());
    const T tolerance = T{4096} * std::numeric_limits<T>::epsilon();

    for (std::size_t i = 0; i < numerator_complex.size(); ++i) {
        if (std::abs(numerator_complex[i].imag()) >
            tolerance * (T{1} + std::abs(numerator_complex[i].real())))
            throw std::runtime_error("numerator is not real");
        output.numerator[i] = zpk.gain * numerator_complex[i].real();
    }
    for (std::size_t i = 0; i < denominator_complex.size(); ++i) {
        if (std::abs(denominator_complex[i].imag()) >
            tolerance * (T{1} + std::abs(denominator_complex[i].real())))
            throw std::runtime_error("denominator is not real");
        output.denominator[i] = denominator_complex[i].real();
    }
    return output;
}

template <Scalar T>
[[nodiscard]] SosDesign<T> second_order_sections(const Zpk<T>& zpk) {
    if (zpk.zeros.size() != zpk.poles.size())
        throw std::invalid_argument("SOS conversion requires equal zero/pole degree");
    auto zero_pairs = detail::pair_roots<T>(zpk.zeros);
    auto pole_pairs = detail::pair_roots<T>(zpk.poles);
    if (zero_pairs.size() != pole_pairs.size())
        throw std::runtime_error("zero/pole section count mismatch");

    SosDesign<T> output;
    output.gain = zpk.gain;
    output.sections.reserve(pole_pairs.size());
    for (std::size_t i = 0; i < pole_pairs.size(); ++i)
        output.sections.push_back(detail::section_from_pairs<T>(zero_pairs[i], pole_pairs[i]));
    return output;
}

}  // namespace signal_processing::iir::design
