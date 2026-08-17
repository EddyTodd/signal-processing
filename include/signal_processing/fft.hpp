#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <numeric>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace signal_processing::fft {

template <typename T>
concept Scalar = std::same_as<T, float> || std::same_as<T, double>;

enum class Direction { forward, inverse };

template <Scalar T>
using Complex = std::complex<T>;

template <Scalar T>
using Vector = std::vector<Complex<T>>;

namespace detail {

[[nodiscard]] constexpr bool is_power_of_two(std::size_t n) noexcept {
    return n != 0 && (n & (n - 1)) == 0;
}

[[nodiscard]] inline std::size_t next_power_of_two(std::size_t n) {
    if (n <= 1) return 1;
    const auto high = std::size_t{1} << (std::numeric_limits<std::size_t>::digits - 1);
    if (n > high) throw std::overflow_error("next power of two overflows size_t");
    --n;
    for (std::size_t shift = 1; shift < std::numeric_limits<std::size_t>::digits; shift <<= 1) {
        n |= n >> shift;
    }
    return n + 1;
}

[[nodiscard]] inline std::size_t ilog2(std::size_t n) {
    if (!is_power_of_two(n)) throw std::invalid_argument("ilog2 requires a power of two");
    std::size_t result = 0;
    while (n > 1) { n >>= 1; ++result; }
    return result;
}

[[nodiscard]] inline bool is_prime(std::size_t n) noexcept {
    if (n < 2) return false;
    if ((n & 1U) == 0U) return n == 2;
    for (std::size_t d = 3; d <= n / d; d += 2) if (n % d == 0) return false;
    return true;
}

[[nodiscard]] inline std::size_t smallest_factor(std::size_t n) noexcept {
    if ((n & 1U) == 0U) return 2;
    for (std::size_t f = 3; f <= n / f; f += 2) if (n % f == 0) return f;
    return n;
}

[[nodiscard]] inline std::size_t mul_mod(std::size_t a, std::size_t b, std::size_t modulus) noexcept {
    std::size_t result = 0;
    a %= modulus;
    while (b != 0) {
        if ((b & 1U) != 0U) result = a >= modulus - result ? a - (modulus - result) : a + result;
        b >>= 1;
        if (b != 0) a = a >= modulus - a ? a - (modulus - a) : a + a;
    }
    return result;
}

[[nodiscard]] inline std::size_t pow_mod(std::size_t base, std::size_t exponent,
                                         std::size_t modulus) noexcept {
    std::size_t result = 1 % modulus;
    while (exponent != 0) {
        if ((exponent & 1U) != 0U) result = mul_mod(result, base, modulus);
        exponent >>= 1;
        if (exponent != 0) base = mul_mod(base, base, modulus);
    }
    return result;
}

[[nodiscard]] inline std::vector<std::size_t> prime_factors(std::size_t n) {
    std::vector<std::size_t> factors;
    for (std::size_t f = 2; f <= n / f; ++f) {
        if (n % f != 0) continue;
        factors.push_back(f);
        while (n % f == 0) n /= f;
    }
    if (n > 1) factors.push_back(n);
    return factors;
}

[[nodiscard]] inline std::size_t primitive_root_prime(std::size_t prime) {
    if (prime == 2) return 1;
    if (!is_prime(prime)) throw std::invalid_argument("primitive root requires a prime length");
    const auto factors = prime_factors(prime - 1);
    for (std::size_t candidate = 2; candidate < prime; ++candidate) {
        bool valid = true;
        for (const auto factor : factors) {
            if (pow_mod(candidate, (prime - 1) / factor, prime) == 1) { valid = false; break; }
        }
        if (valid) return candidate;
    }
    throw std::runtime_error("failed to find primitive root");
}

[[nodiscard]] inline std::size_t modular_inverse(std::size_t a, std::size_t modulus) {
    if (modulus == 0) throw std::invalid_argument("modular inverse modulus must be nonzero");
    using Signed = std::int64_t;
    if (a > static_cast<std::size_t>(std::numeric_limits<Signed>::max()) ||
        modulus > static_cast<std::size_t>(std::numeric_limits<Signed>::max())) {
        throw std::overflow_error("modular inverse input exceeds supported signed range");
    }
    Signed t = 0, next_t = 1;
    Signed r = static_cast<Signed>(modulus), next_r = static_cast<Signed>(a % modulus);
    while (next_r != 0) {
        const Signed q = r / next_r;
        std::tie(t, next_t) = std::pair<Signed, Signed>{next_t, t - q * next_t};
        std::tie(r, next_r) = std::pair<Signed, Signed>{next_r, r - q * next_r};
    }
    if (r != 1) throw std::invalid_argument("modular inverse does not exist");
    if (t < 0) t += static_cast<Signed>(modulus);
    return static_cast<std::size_t>(t);
}

[[nodiscard]] inline std::pair<std::size_t, std::size_t> coprime_factor_split(std::size_t n) noexcept {
    std::pair<std::size_t, std::size_t> best{0, 0};
    std::size_t best_gap = n;
    for (std::size_t a = 2; a <= n / a; ++a) {
        if (n % a != 0) continue;
        const auto b = n / a;
        if (std::gcd(a, b) != 1) continue;
        const auto gap = b - a;
        if (gap < best_gap) { best = {a, b}; best_gap = gap; }
    }
    return best;
}

template <Scalar T>
[[nodiscard]] inline Complex<T> root(T angle) {
    return {std::cos(angle), std::sin(angle)};
}

template <Scalar T>
[[nodiscard]] inline T sign(Direction direction) noexcept {
    return direction == Direction::forward ? T{-1} : T{1};
}

template <Scalar T>
inline void normalize_inverse(std::span<Complex<T>> data, Direction direction) {
    if (direction != Direction::inverse || data.empty()) return;
    const T scale = T{1} / static_cast<T>(data.size());
    for (auto& value : data) value *= scale;
}

template <Scalar T>
[[nodiscard]] inline T chirp_phase_ratio(std::size_t k, std::size_t n) {
    if (n > std::numeric_limits<std::size_t>::max() / 2)
        throw std::length_error("chirp phase modulus overflow");
    const std::size_t period = 2 * n;
    const std::size_t residue = mul_mod(k % period, k % period, period);
    return static_cast<T>(residue) / static_cast<T>(n);
}

template <Scalar T>
[[nodiscard]] Vector<T> dft_unnormalized(std::span<const Complex<T>> input, Direction direction) {
    const auto n = input.size();
    Vector<T> output(n);
    if (n == 0) return output;
    const T tau = T{2} * std::numbers::pi_v<T>;
    const T s = sign<T>(direction);
    for (std::size_t k = 0; k < n; ++k) {
        Complex<T> sum{};
        for (std::size_t t = 0; t < n; ++t) {
            const T angle = s * tau * static_cast<T>(k) * static_cast<T>(t) / static_cast<T>(n);
            sum += input[t] * root<T>(angle);
        }
        output[k] = sum;
    }
    return output;
}

template <Scalar T>
void recursive_core(Vector<T>& data, Direction direction) {
    const auto n = data.size();
    if (n <= 1) return;
    Vector<T> even(n / 2), odd(n / 2);
    for (std::size_t i = 0; i < n / 2; ++i) { even[i] = data[2 * i]; odd[i] = data[2 * i + 1]; }
    recursive_core(even, direction);
    recursive_core(odd, direction);
    Complex<T> w{1, 0};
    const auto step = root<T>(sign<T>(direction) * T{2} * std::numbers::pi_v<T> / static_cast<T>(n));
    for (std::size_t k = 0; k < n / 2; ++k) {
        const auto v = w * odd[k];
        data[k] = even[k] + v;
        data[k + n / 2] = even[k] - v;
        w *= step;
    }
}

template <Scalar T>
void radix4_core(Vector<T>& data, Direction direction) {
    const auto n = data.size();
    if (n <= 1) return;
    if (n == 2) {
        const auto a = data[0], b = data[1];
        data[0] = a + b; data[1] = a - b; return;
    }
    const auto quarter = n / 4;
    std::array<Vector<T>, 4> sub{Vector<T>(quarter), Vector<T>(quarter), Vector<T>(quarter), Vector<T>(quarter)};
    for (std::size_t r = 0; r < 4; ++r) for (std::size_t j = 0; j < quarter; ++j) sub[r][j] = data[4 * j + r];
    for (auto& s : sub) radix4_core(s, direction);
    const T sg = sign<T>(direction);
    const Complex<T> imag_unit{0, sg};
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k = 0; k < quarter; ++k) {
        const auto b = sub[1][k] * root<T>(sg * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto c = sub[2][k] * root<T>(sg * T{2} * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto d = sub[3][k] * root<T>(sg * T{3} * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto even_sum = sub[0][k] + c;
        const auto even_diff = sub[0][k] - c;
        const auto odd_sum = b + d;
        const auto odd_diff = (b - d) * imag_unit;
        data[k] = even_sum + odd_sum;
        data[k + quarter] = even_diff + odd_diff;
        data[k + 2 * quarter] = even_sum - odd_sum;
        data[k + 3 * quarter] = even_diff - odd_diff;
    }
}

template <Scalar T>
void split_core(Vector<T>& data, Direction direction) {
    const auto n = data.size();
    if (n <= 1) return;
    if (n == 2) {
        const auto a = data[0], b = data[1]; data[0] = a + b; data[1] = a - b; return;
    }
    const auto quarter = n / 4;
    Vector<T> even(n / 2), odd1(quarter), odd3(quarter);
    for (std::size_t j = 0; j < n / 2; ++j) even[j] = data[2 * j];
    for (std::size_t j = 0; j < quarter; ++j) { odd1[j] = data[4 * j + 1]; odd3[j] = data[4 * j + 3]; }
    split_core(even, direction); split_core(odd1, direction); split_core(odd3, direction);
    const T sg = sign<T>(direction);
    const Complex<T> imag_unit{0, sg};
    const T tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k = 0; k < quarter; ++k) {
        const auto a = odd1[k] * root<T>(sg * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto b = odd3[k] * root<T>(sg * T{3} * tau * static_cast<T>(k) / static_cast<T>(n));
        const auto sum = a + b, diff = a - b;
        data[k] = even[k] + sum;
        data[k + n / 2] = even[k] - sum;
        data[k + quarter] = even[k + quarter] + imag_unit * diff;
        data[k + 3 * quarter] = even[k + quarter] - imag_unit * diff;
    }
}

template <Scalar T>
[[nodiscard]] T modified_split_scale(std::size_t n, std::size_t k) {
    if (n <= 4) return T{1};
    const auto q = n / 4;
    const auto rk = k % q;
    const auto sub = modified_split_scale<T>(q, rk);
    const auto angle = T{2} * std::numbers::pi_v<T> * static_cast<T>(rk) / static_cast<T>(n);
    return sub * (rk <= n / 8 ? std::cos(angle) : std::sin(angle));
}

template <Scalar T>
[[nodiscard]] Complex<T> modified_twiddle_mul(const Complex<T>& value, std::size_t n,
                                               std::size_t k, Direction direction,
                                               bool conjugate) {
    const auto rk = k % (n / 4);
    const auto theta = T{2} * std::numbers::pi_v<T> * static_cast<T>(rk) / static_cast<T>(n);
    T is = direction == Direction::forward ? T{-1} : T{1};
    if (conjugate) is = -is;
    const T re = value.real(), im = value.imag();
    if (rk <= n / 8) {
        const T tangent = std::tan(theta);
        return {re - is * im * tangent, im + is * re * tangent};
    }
    const T cotangent = T{1} / std::tan(theta);
    return {re * cotangent - is * im, im * cotangent + is * re};
}

template <Scalar T>
void modified_split_scaled_core(const Vector<T>& input, Vector<T>& scaled, Direction direction) {
    const auto n = input.size();
    scaled.resize(n);
    if (n <= 4) { scaled = dft_unnormalized<T>(input, direction); return; }
    Vector<T> even_in(n / 2), odd1_in(n / 4), odd3_in(n / 4);
    for (std::size_t j = 0; j < n / 2; ++j) even_in[j] = input[2 * j];
    for (std::size_t j = 0; j < n / 4; ++j) {
        odd1_in[j] = input[4 * j + 1];
        odd3_in[j] = input[(4 * j + n - 1) % n];
    }
    Vector<T> even, odd1, odd3;
    modified_split_scaled_core(even_in, even, direction);
    modified_split_scaled_core(odd1_in, odd1, direction);
    modified_split_scaled_core(odd3_in, odd3, direction);
    const auto q = n / 4;
    for (std::size_t k = 0; k < q; ++k) {
        const T scale = modified_split_scale<T>(n, k);
        const auto first = modified_twiddle_mul<T>(odd1[k], n, k, direction, false);
        const auto third = modified_twiddle_mul<T>(odd3[k], n, k, direction, true);
        const auto sum = first + third, diff = first - third;
        const auto e0 = even[k] * (modified_split_scale<T>(n / 2, k) / scale);
        const auto e1 = even[k + q] * (modified_split_scale<T>(n / 2, k + q) / scale);
        const Complex<T> minus_i = direction == Direction::forward ? Complex<T>{0, -1} : Complex<T>{0, 1};
        scaled[k] = e0 + sum;
        scaled[k + q] = e1 + minus_i * diff;
        scaled[k + 2 * q] = e0 - sum;
        scaled[k + 3 * q] = e1 - minus_i * diff;
    }
}

template <Scalar T>
[[nodiscard]] Vector<T> mixed_core(std::span<const Complex<T>> input, Direction direction) {
    const auto n = input.size();
    if (n <= 1) return Vector<T>(input.begin(), input.end());
    const auto radix = smallest_factor(n);
    if (radix == n) return dft_unnormalized<T>(input, direction);
    const auto m = n / radix;
    std::vector<Vector<T>> sub(radix, Vector<T>(m));
    for (std::size_t q = 0; q < radix; ++q) {
        for (std::size_t j = 0; j < m; ++j) sub[q][j] = input[radix * j + q];
        sub[q] = mixed_core<T>(sub[q], direction);
    }
    Vector<T> output(n);
    const T sg = sign<T>(direction), tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t k0 = 0; k0 < m; ++k0) {
        for (std::size_t k1 = 0; k1 < radix; ++k1) {
            const auto k = k0 + m * k1;
            Complex<T> sum{};
            for (std::size_t q = 0; q < radix; ++q) {
                const T angle = sg * tau * static_cast<T>(k) * static_cast<T>(q) / static_cast<T>(n);
                sum += sub[q][k0] * root<T>(angle);
            }
            output[k] = sum;
        }
    }
    return output;
}

template <Scalar T>
[[nodiscard]] Vector<T> circular_convolution_fft(Vector<T> a, Vector<T> b);

}  // namespace detail

[[nodiscard]] constexpr bool is_power_of_two(std::size_t n) noexcept { return detail::is_power_of_two(n); }

template <Scalar T>
[[nodiscard]] inline Vector<T> dft(std::span<const Complex<T>> input,
                                   Direction direction = Direction::forward) {
    auto output = detail::dft_unnormalized<T>(input, direction);
    detail::normalize_inverse<T>(output, direction);
    return output;
}

template <Scalar T>
inline void radix2_inplace(std::span<Complex<T>> data, Direction direction = Direction::forward) {
    const auto n = data.size();
    if (n <= 1) return;
    if (!detail::is_power_of_two(n)) throw std::invalid_argument("radix2_inplace requires a power-of-two length");
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(data[i], data[j]);
    }
    const T sg = detail::sign<T>(direction), tau = T{2} * std::numbers::pi_v<T>;
    for (std::size_t length = 2; length <= n;) {
        const auto step = detail::root<T>(sg * tau / static_cast<T>(length));
        for (std::size_t base = 0; base < n; base += length) {
            Complex<T> w{1, 0};
            for (std::size_t off = 0; off < length / 2; ++off) {
                const auto u = data[base + off];
                const auto v = data[base + off + length / 2] * w;
                data[base + off] = u + v;
                data[base + off + length / 2] = u - v;
                w *= step;
            }
        }
        if (length == n) break;
        length <<= 1;
    }
    detail::normalize_inverse<T>(data, direction);
}

template <Scalar T>
[[nodiscard]] inline Vector<T> radix2(std::span<const Complex<T>> input,
                                      Direction direction = Direction::forward) {
    Vector<T> output(input.begin(), input.end());
    radix2_inplace<T>(output, direction);
    return output;
}

template <Scalar T>
[[nodiscard]] inline Vector<T> radix2_recursive(std::span<const Complex<T>> input,
                                                Direction direction = Direction::forward) {
    if (!input.empty() && !detail::is_power_of_two(input.size()))
        throw std::invalid_argument("radix2_recursive requires a power-of-two length");
    Vector<T> output(input.begin(), input.end());
    detail::recursive_core(output, direction);
    detail::normalize_inverse<T>(output, direction);
    return output;
}

template <Scalar T>
[[nodiscard]] inline Vector<T> stockham(std::span<const Complex<T>> input,
                                        Direction direction = Direction::forward) {
    const auto n = input.size();
    if (n <= 1) return Vector<T>(input.begin(), input.end());
    if (!detail::is_power_of_two(n)) throw std::invalid_argument("stockham requires a power-of-two length");
    Vector<T> source(input.begin(), input.end()), destination(n);
    const T sg = detail::sign<T>(direction);
    for (std::size_t m = 1; m < n; m <<= 1) {
        for (std::size_t j = 0; j < n / 2; ++j) {
            const auto k = j & (m - 1);
            const auto w = detail::root<T>(sg * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(m));
            const auto lower = source[j + n / 2] * w;
            const auto out = 2 * j - k;
            destination[out] = source[j] + lower;
            destination[out + m] = source[j] - lower;
        }
        source.swap(destination);
    }
    detail::normalize_inverse<T>(source, direction);
    return source;
}

template <Scalar T>
[[nodiscard]] inline Vector<T> radix4(std::span<const Complex<T>> input,
                                      Direction direction = Direction::forward) {
    if (!input.empty() && !detail::is_power_of_two(input.size()))
        throw std::invalid_argument("radix4 requires a power-of-two length");
    Vector<T> output(input.begin(), input.end());
    detail::radix4_core(output, direction);
    detail::normalize_inverse<T>(output, direction);
    return output;
}

template <Scalar T>
[[nodiscard]] inline Vector<T> split_radix(std::span<const Complex<T>> input,
                                           Direction direction = Direction::forward) {
    if (!input.empty() && !detail::is_power_of_two(input.size()))
        throw std::invalid_argument("split_radix requires a power-of-two length");
    Vector<T> output(input.begin(), input.end());
    detail::split_core(output, direction);
    detail::normalize_inverse<T>(output, direction);
    return output;
}

template <Scalar T>
[[nodiscard]] inline Vector<T> modified_split_radix(std::span<const Complex<T>> input,
                                                    Direction direction = Direction::forward) {
    if (!input.empty() && !detail::is_power_of_two(input.size()))
        throw std::invalid_argument("modified_split_radix requires a power-of-two length");
    if (input.size() <= 1) return Vector<T>(input.begin(), input.end());
    Vector<T> in(input.begin(), input.end()), scaled;
    detail::modified_split_scaled_core(in, scaled, direction);
    for (std::size_t k = 0; k < scaled.size(); ++k)
        scaled[k] *= detail::modified_split_scale<T>(scaled.size(), k);
    detail::normalize_inverse<T>(scaled, direction);
    return scaled;
}

template <Scalar T>
[[nodiscard]] inline Vector<T> mixed_radix(std::span<const Complex<T>> input,
                                           Direction direction = Direction::forward) {
    auto output = detail::mixed_core<T>(input, direction);
    detail::normalize_inverse<T>(output, direction);
    return output;
}

template <Scalar T>
[[nodiscard]] inline Vector<T> good_thomas(std::span<const Complex<T>> input,
                                           std::size_t factor_a = 0, std::size_t factor_b = 0,
                                           Direction direction = Direction::forward) {
    const auto n = input.size();
    if (n <= 1) return Vector<T>(input.begin(), input.end());
    if (factor_a == 0 || factor_b == 0) std::tie(factor_a, factor_b) = detail::coprime_factor_split(n);
    if (factor_a <= 1 || factor_b <= 1 || factor_a * factor_b != n || std::gcd(factor_a, factor_b) != 1)
        throw std::invalid_argument("good_thomas requires a coprime factorization N=a*b");
    const auto inv_b_a = detail::modular_inverse(factor_b % factor_a, factor_a);
    const auto inv_a_b = detail::modular_inverse(factor_a % factor_b, factor_b);
    Vector<T> matrix(n);
    for (std::size_t n1 = 0; n1 < factor_a; ++n1) {
        for (std::size_t n2 = 0; n2 < factor_b; ++n2) {
            const auto index = (detail::mul_mod(detail::mul_mod(n1, factor_b, n), inv_b_a, n) +
                                detail::mul_mod(detail::mul_mod(n2, factor_a, n), inv_a_b, n)) % n;
            matrix[n1 * factor_b + n2] = input[index];
        }
    }
    for (std::size_t n1 = 0; n1 < factor_a; ++n1) {
        auto row = mixed_radix<T>(std::span<const Complex<T>>(matrix.data() + n1 * factor_b, factor_b), direction);
        std::copy(row.begin(), row.end(), matrix.begin() + static_cast<std::ptrdiff_t>(n1 * factor_b));
    }
    for (std::size_t k2 = 0; k2 < factor_b; ++k2) {
        Vector<T> column(factor_a);
        for (std::size_t n1 = 0; n1 < factor_a; ++n1) column[n1] = matrix[n1 * factor_b + k2];
        column = mixed_radix<T>(column, direction);
        for (std::size_t k1 = 0; k1 < factor_a; ++k1) matrix[k1 * factor_b + k2] = column[k1];
    }
    Vector<T> output(n);
    for (std::size_t k1 = 0; k1 < factor_a; ++k1)
        for (std::size_t k2 = 0; k2 < factor_b; ++k2)
            output[(k1 * factor_b + k2 * factor_a) % n] = matrix[k1 * factor_b + k2];
    return output;
}

template <Scalar T>
[[nodiscard]] inline Vector<T> bluestein(std::span<const Complex<T>> input,
                                         Direction direction = Direction::forward) {
    const auto n = input.size();
    if (n <= 1) return Vector<T>(input.begin(), input.end());
    if (n > (std::numeric_limits<std::size_t>::max() / 2) + 1)
        throw std::length_error("bluestein workspace overflow");
    const auto m = detail::next_power_of_two(2 * n - 1);
    Vector<T> a(m), b(m);
    const T sg = detail::sign<T>(direction);
    for (std::size_t k = 0; k < n; ++k) {
        const T phase = detail::chirp_phase_ratio<T>(k, n);
        const auto chirp = detail::root<T>(sg * std::numbers::pi_v<T> * phase);
        a[k] = input[k] * chirp;
        b[k] = std::conj(chirp);
        if (k != 0) b[m - k] = std::conj(chirp);
    }
    radix2_inplace<T>(a, Direction::forward);
    radix2_inplace<T>(b, Direction::forward);
    for (std::size_t i = 0; i < m; ++i) a[i] *= b[i];
    radix2_inplace<T>(a, Direction::inverse);
    Vector<T> output(n);
    for (std::size_t k = 0; k < n; ++k) {
        const T phase = detail::chirp_phase_ratio<T>(k, n);
        output[k] = a[k] * detail::root<T>(sg * std::numbers::pi_v<T> * phase);
    }
    detail::normalize_inverse<T>(output, direction);
    return output;
}

namespace detail {
template <Scalar T>
[[nodiscard]] inline Vector<T> circular_convolution_fft(Vector<T> a, Vector<T> b) {
    if (a.size() != b.size()) throw std::invalid_argument("circular convolution size mismatch");
    const auto n = a.size();
    if (n == 0) return {};
    const auto m = next_power_of_two(2 * n - 1);
    a.resize(m); b.resize(m);
    radix2_inplace<T>(a, Direction::forward); radix2_inplace<T>(b, Direction::forward);
    for (std::size_t i = 0; i < m; ++i) a[i] *= b[i];
    radix2_inplace<T>(a, Direction::inverse);
    Vector<T> output(n);
    for (std::size_t i = 0; i < 2 * n - 1; ++i) output[i % n] += a[i];
    return output;
}
} // namespace detail

template <Scalar T>
[[nodiscard]] inline Vector<T> rader(std::span<const Complex<T>> input,
                                     Direction direction = Direction::forward) {
    const auto n = input.size();
    if (n <= 2) return dft<T>(input, direction);
    if (!detail::is_prime(n)) throw std::invalid_argument("rader requires a prime length");
    const auto generator = detail::primitive_root_prime(n);
    const auto l = n - 1;
    Vector<T> a(l), b(l);
    std::vector<std::size_t> powers(l);
    powers[0] = 1;
    for (std::size_t q = 1; q < l; ++q) powers[q] = detail::mul_mod(powers[q - 1], generator, n);
    const T sg = detail::sign<T>(direction);
    for (std::size_t q = 0; q < l; ++q) {
        a[q] = input[powers[(l - q) % l]];
        const T angle = sg * T{2} * std::numbers::pi_v<T> * static_cast<T>(powers[q]) / static_cast<T>(n);
        b[q] = detail::root<T>(angle);
    }
    const auto convolution = detail::circular_convolution_fft(std::move(a), std::move(b));
    Vector<T> output(n);
    output[0] = std::accumulate(input.begin(), input.end(), Complex<T>{});
    for (std::size_t q = 0; q < l; ++q) output[powers[q]] = input[0] + convolution[q];
    detail::normalize_inverse<T>(output, direction);
    return output;
}

}  // namespace signal_processing::fft

#include "signal_processing/detail/fft_plans.hpp"
#include "signal_processing/detail/fft_simd.hpp"
