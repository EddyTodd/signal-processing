#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <memory>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace signal_processing::fft {

[[nodiscard]] inline constexpr bool small_codelet_supported(std::size_t radix) noexcept {
    return radix == 2 || radix == 3 || radix == 4 || radix == 5 || radix == 7;
}

[[nodiscard]] inline std::size_t choose_small_radix(std::size_t n) noexcept {
    for (const auto r : {4U, 2U, 3U, 5U, 7U}) if (n % r == 0) return r;
    return 0;
}

template <Scalar T = double>
class SmallDftCodelet {
public:
    explicit SmallDftCodelet(std::size_t radix) : radix_(radix) {
        if (!small_codelet_supported(radix_))
            throw std::invalid_argument("SmallDftCodelet supports radix 2, 3, 4, 5, or 7");
        if (radix_ == 3 || radix_ == 5 || radix_ == 7) {
            roots_.resize(radix_ * radix_);
            const T tau = T{2} * std::numbers::pi_v<T>;
            for (std::size_t k = 0; k < radix_; ++k)
                for (std::size_t q = 0; q < radix_; ++q)
                    roots_[k * radix_ + q] = detail::root<T>(
                        -tau * static_cast<T>(k * q) / static_cast<T>(radix_));
        }
    }

    [[nodiscard]] std::size_t radix() const noexcept { return radix_; }
    [[nodiscard]] std::size_t stored_roots() const noexcept { return roots_.size(); }

    void execute(std::span<Complex<T>> values,
                 Direction direction = Direction::forward) const {
        if (values.size() != radix_)
            throw std::invalid_argument("SmallDftCodelet buffer size mismatch");
        if (radix_ == 2) {
            const auto a = values[0], b = values[1];
            values[0] = a + b; values[1] = a - b;
        } else if (radix_ == 4) {
            const auto a = values[0], b = values[1], c = values[2], d = values[3];
            const auto e0 = a + c, e1 = a - c, o0 = b + d, diff = b - d;
            const Complex<T> i_term = direction == Direction::forward
                ? Complex<T>{diff.imag(), -diff.real()}
                : Complex<T>{-diff.imag(), diff.real()};
            values[0] = e0 + o0; values[2] = e0 - o0;
            values[1] = e1 + i_term; values[3] = e1 - i_term;
        } else {
            std::array<Complex<T>, 7> input{};
            for (std::size_t i = 0; i < radix_; ++i) input[i] = values[i];
            for (std::size_t k = 0; k < radix_; ++k) {
                Complex<T> sum{};
                for (std::size_t q = 0; q < radix_; ++q) {
                    auto w = roots_[k * radix_ + q];
                    if (direction == Direction::inverse) w = std::conj(w);
                    sum += input[q] * w;
                }
                values[k] = sum;
            }
        }
        if (direction == Direction::inverse) {
            const T scale = T{1} / static_cast<T>(radix_);
            for (auto& value : values) value *= scale;
        }
    }

private:
    std::size_t radix_{};
    Vector<T> roots_;
};
SmallDftCodelet(std::size_t) -> SmallDftCodelet<double>;

template <Scalar T = double>
class Radix2Plan {
public:
    explicit Radix2Plan(std::size_t n)
        : n_(checked_length(n)), bit_reverse_(n_), twiddles_(n_ / 2) {
        const auto bits = detail::ilog2(n_);
        for (std::size_t i = 0; i < n_; ++i) {
            auto x = i; std::size_t reversed = 0;
            for (std::size_t b = 0; b < bits; ++b) { reversed = (reversed << 1) | (x & 1U); x >>= 1; }
            bit_reverse_[i] = reversed;
        }
        for (std::size_t k = 0; k < n_ / 2; ++k) {
            const T angle = -T{2} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n_);
            twiddles_[k] = detail::root<T>(angle);
        }
    }

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return 0; }
    [[nodiscard]] std::size_t stored_twiddles() const noexcept { return twiddles_.size(); }
    [[nodiscard]] std::size_t stored_indices() const noexcept { return bit_reverse_.size(); }

    void forward_inplace(std::span<Complex<T>> data) const { execute(data, Direction::forward); }
    void inverse_inplace(std::span<Complex<T>> data) const { execute(data, Direction::inverse); }
    void forward(std::span<const Complex<T>> input, std::span<Complex<T>> output) const {
        copy_checked(input, output); forward_inplace(output);
    }
    void inverse(std::span<const Complex<T>> input, std::span<Complex<T>> output) const {
        copy_checked(input, output); inverse_inplace(output);
    }

private:
    [[nodiscard]] static std::size_t checked_length(std::size_t n) {
        if (n == 0 || !detail::is_power_of_two(n))
            throw std::invalid_argument("Radix2Plan requires a power-of-two length >= 1");
        return n;
    }
    void copy_checked(std::span<const Complex<T>> input, std::span<Complex<T>> output) const {
        if (input.size() != n_ || output.size() != n_)
            throw std::invalid_argument("Radix2Plan buffer size mismatch");
        std::copy(input.begin(), input.end(), output.begin());
    }
    void execute(std::span<Complex<T>> data, Direction direction) const {
        if (data.size() != n_) throw std::invalid_argument("Radix2Plan buffer size mismatch");
        if (n_ <= 1) return;
        for (std::size_t i = 0; i < n_; ++i)
            if (i < bit_reverse_[i]) std::swap(data[i], data[bit_reverse_[i]]);
        for (std::size_t len = 2; len <= n_;) {
            const auto stride = n_ / len;
            for (std::size_t base = 0; base < n_; base += len) {
                for (std::size_t j = 0; j < len / 2; ++j) {
                    auto w = twiddles_[j * stride];
                    if (direction == Direction::inverse) w = std::conj(w);
                    const auto u = data[base + j];
                    const auto v = data[base + j + len / 2] * w;
                    data[base + j] = u + v;
                    data[base + j + len / 2] = u - v;
                }
            }
            if (len == n_) break;
            len <<= 1;
        }
        detail::normalize_inverse<T>(data, direction);
    }

    std::size_t n_{};
    std::vector<std::size_t> bit_reverse_;
    Vector<T> twiddles_;
};
Radix2Plan(std::size_t) -> Radix2Plan<double>;

template <Scalar T = double>
class RealRadix2Plan {
public:
    explicit RealRadix2Plan(std::size_t n)
        : n_(checked_length(n)), half_(n_ / 2), half_plan_(n_ == 1 ? 1 : n_ / 2),
          post_twiddles_(n_ == 1 ? 1 : n_ / 2 + 1) {
        if (n_ == 1) { post_twiddles_[0] = {1, 0}; return; }
        for (std::size_t k = 0; k <= half_; ++k) {
            const T angle = -T{2} * std::numbers::pi_v<T> * static_cast<T>(k) / static_cast<T>(n_);
            post_twiddles_[k] = detail::root<T>(angle);
        }
    }

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t spectrum_size() const noexcept { return n_ == 1 ? 1 : half_ + 1; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return n_ == 1 ? 0 : half_; }

    void forward(std::span<const T> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const {
        check(input.size(), output.size(), scratch.size());
        if (n_ == 1) { output[0] = {input[0], 0}; return; }
        for (std::size_t j = 0; j < half_; ++j) scratch[j] = {input[2 * j], input[2 * j + 1]};
        half_plan_.forward_inplace(scratch.first(half_));
        for (std::size_t k = 0; k <= half_; ++k) {
            const auto a = scratch[k % half_];
            const auto b = std::conj(scratch[(half_ - k) % half_]);
            const auto even = T{0.5} * (a + b);
            const auto odd = Complex<T>{0, T{-0.5}} * (a - b);
            output[k] = even + post_twiddles_[k] * odd;
        }
    }

    void inverse(std::span<const Complex<T>> input, std::span<T> output,
                 std::span<Complex<T>> scratch) const {
        check(output.size(), input.size(), scratch.size());
        if (n_ == 1) { output[0] = input[0].real(); return; }
        for (std::size_t k = 0; k < half_; ++k) {
            const auto xk = input[k];
            const auto mirror = std::conj(input[k != 0 ? half_ - k : half_]);
            const auto even = T{0.5} * (xk + mirror);
            const auto odd = T{0.5} * (xk - mirror) / post_twiddles_[k];
            scratch[k] = even + Complex<T>{0, 1} * odd;
        }
        half_plan_.inverse_inplace(scratch.first(half_));
        for (std::size_t j = 0; j < half_; ++j) {
            output[2 * j] = scratch[j].real();
            output[2 * j + 1] = scratch[j].imag();
        }
    }

private:
    [[nodiscard]] static std::size_t checked_length(std::size_t n) {
        if (n == 0 || !detail::is_power_of_two(n))
            throw std::invalid_argument("RealRadix2Plan requires a power-of-two length >= 1");
        return n;
    }
    void check(std::size_t time_size, std::size_t freq_size, std::size_t scratch_size) const {
        if (time_size != n_ || freq_size != spectrum_size() || scratch_size < this->scratch_size())
            throw std::invalid_argument("RealRadix2Plan buffer size mismatch");
    }
    std::size_t n_{}, half_{};
    Radix2Plan<T> half_plan_;
    Vector<T> post_twiddles_;
};
RealRadix2Plan(std::size_t) -> RealRadix2Plan<double>;

template <Scalar T = double>
class MixedRadixPlan {
    struct Node {
        std::size_t n{}, radix{}, m{};
        bool leaf{};
        Vector<T> twiddles;
        std::unique_ptr<SmallDftCodelet<T>> codelet;
        std::unique_ptr<Node> child;
    };

public:
    explicit MixedRadixPlan(std::size_t n) : n_(n), root_(build(n)) {
        if (n == 0) throw std::invalid_argument("MixedRadixPlan requires length >= 1");
    }
    MixedRadixPlan(MixedRadixPlan&&) noexcept = default;
    MixedRadixPlan& operator=(MixedRadixPlan&&) noexcept = default;
    MixedRadixPlan(const MixedRadixPlan&) = delete;
    MixedRadixPlan& operator=(const MixedRadixPlan&) = delete;

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return n_; }
    [[nodiscard]] std::size_t top_radix() const noexcept { return root_->radix; }

    void forward_inplace(std::span<Complex<T>> data, std::span<Complex<T>> scratch) const {
        execute(data, scratch, Direction::forward);
    }
    void inverse_inplace(std::span<Complex<T>> data, std::span<Complex<T>> scratch) const {
        execute(data, scratch, Direction::inverse);
    }
    void forward(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const {
        copy_checked(input, output, scratch); forward_inplace(output, scratch);
    }
    void inverse(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const {
        copy_checked(input, output, scratch); inverse_inplace(output, scratch);
    }

private:
    [[nodiscard]] static std::unique_ptr<Node> build(std::size_t n) {
        if (n == 0) throw std::invalid_argument("MixedRadixPlan requires length >= 1");
        auto node = std::make_unique<Node>(); node->n = n;
        if (n == 1) { node->leaf = true; node->radix = 1; node->m = 1; return node; }
        if (small_codelet_supported(n)) {
            node->leaf = true; node->radix = n; node->m = 1;
            node->codelet = std::make_unique<SmallDftCodelet<T>>(n); return node;
        }
        const auto radix = choose_small_radix(n);
        if (radix == 0) {
            node->leaf = true; node->radix = n; node->m = 1;
            if (n > std::numeric_limits<std::size_t>::max() / n)
                throw std::length_error("MixedRadixPlan direct leaf table overflows size_t");
            node->twiddles.resize(n * n);
            for (std::size_t k = 0; k < n; ++k) {
                for (std::size_t q = 0; q < n; ++q) {
                    const std::size_t phase = detail::mul_mod(k, q, n);
                    node->twiddles[k * n + q] = detail::root<T>(
                        -T{2} * std::numbers::pi_v<T> * static_cast<T>(phase) /
                        static_cast<T>(n));
                }
            }
            return node;
        }
        node->radix = radix; node->m = n / radix;
        node->codelet = std::make_unique<SmallDftCodelet<T>>(radix);
        node->child = build(node->m);
        node->twiddles.resize(node->m * radix);
        for (std::size_t k0 = 0; k0 < node->m; ++k0)
            for (std::size_t q = 0; q < radix; ++q)
                node->twiddles[k0 * radix + q] = detail::root<T>(
                    -T{2} * std::numbers::pi_v<T> * static_cast<T>(k0 * q) / static_cast<T>(n));
        return node;
    }

    static void direct_leaf(const Node& node, std::span<Complex<T>> data,
                            std::span<Complex<T>> scratch, Direction direction) {
        const auto n = node.n;
        for (std::size_t k = 0; k < n; ++k) {
            Complex<T> sum{};
            for (std::size_t q = 0; q < n; ++q) {
                auto w = node.twiddles[k * n + q];
                if (direction == Direction::inverse) w = std::conj(w);
                sum += data[q] * w;
            }
            scratch[k] = sum;
        }
        if (direction == Direction::inverse) {
            const T scale = T{1} / static_cast<T>(n);
            for (std::size_t k = 0; k < n; ++k) scratch[k] *= scale;
        }
        const auto result = scratch.first(n);
        std::copy(result.begin(), result.end(), data.begin());
    }

    static void execute_node(const Node& node, std::span<Complex<T>> data,
                             std::span<Complex<T>> scratch, Direction direction) {
        if (node.n == 1) return;
        if (node.leaf) {
            if (node.codelet) node.codelet->execute(data.first(node.n), direction);
            else direct_leaf(node, data.first(node.n), scratch.first(node.n), direction);
            return;
        }
        const auto r = node.radix, m = node.m;
        for (std::size_t q = 0; q < r; ++q)
            for (std::size_t j = 0; j < m; ++j)
                scratch[q * m + j] = data[r * j + q];
        for (std::size_t q = 0; q < r; ++q)
            execute_node(*node.child, scratch.subspan(q * m, m), data.subspan(q * m, m), direction);
        std::array<Complex<T>, 7> values{};
        for (std::size_t k0 = 0; k0 < m; ++k0) {
            for (std::size_t q = 0; q < r; ++q) {
                auto w = node.twiddles[k0 * r + q];
                if (direction == Direction::inverse) w = std::conj(w);
                values[q] = scratch[q * m + k0] * w;
            }
            node.codelet->execute(std::span<Complex<T>>(values.data(), r), direction);
            for (std::size_t k1 = 0; k1 < r; ++k1) data[k0 + m * k1] = values[k1];
        }
    }

    void execute(std::span<Complex<T>> data, std::span<Complex<T>> scratch,
                 Direction direction) const {
        if (data.size() != n_ || scratch.size() < n_)
            throw std::invalid_argument("MixedRadixPlan buffer size mismatch");
        execute_node(*root_, data, scratch.first(n_), direction);
    }
    void copy_checked(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                      std::span<Complex<T>> scratch) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < n_)
            throw std::invalid_argument("MixedRadixPlan buffer size mismatch");
        std::copy(input.begin(), input.end(), output.begin());
    }

    std::size_t n_{};
    std::unique_ptr<Node> root_;
};
MixedRadixPlan(std::size_t) -> MixedRadixPlan<double>;

template <Scalar T = double>
class GoodThomasPlan {
    struct ValidatedFactors {
        std::size_t a{};
        std::size_t b{};
        std::size_t scratch{};
    };

public:
    explicit GoodThomasPlan(std::size_t n)
        : GoodThomasPlan(n, validate(n, detail::coprime_factor_split(n))) {}
    GoodThomasPlan(std::size_t n, std::size_t a, std::size_t b)
        : GoodThomasPlan(n, validate(n, {a, b})) {}

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t factor_a() const noexcept { return a_; }
    [[nodiscard]] std::size_t factor_b() const noexcept { return b_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return scratch_size_; }
    [[nodiscard]] std::size_t twiddle_count() const noexcept { return 0; }

    void forward_inplace(std::span<Complex<T>> data, std::span<Complex<T>> scratch) const {
        execute(data, scratch, Direction::forward);
    }
    void inverse_inplace(std::span<Complex<T>> data, std::span<Complex<T>> scratch) const {
        execute(data, scratch, Direction::inverse);
    }
    void forward(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const {
        copy_checked(input, output, scratch); forward_inplace(output, scratch);
    }
    void inverse(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const {
        copy_checked(input, output, scratch); inverse_inplace(output, scratch);
    }

private:
    [[nodiscard]] static ValidatedFactors validate(
        std::size_t n, std::pair<std::size_t, std::size_t> factors) {
        const auto [a, b] = factors;
        if (a <= 1 || b <= 1 || n % a != 0 || n / a != b || std::gcd(a, b) != 1)
            throw std::invalid_argument("GoodThomasPlan requires coprime N=a*b factors");
        const std::size_t maximum = std::max(a, b);
        if (maximum > (std::numeric_limits<std::size_t>::max() - n) / 2)
            throw std::length_error("GoodThomasPlan scratch size overflows size_t");
        return {a, b, n + 2 * maximum};
    }

    GoodThomasPlan(std::size_t n, ValidatedFactors factors)
        : n_(n), a_(factors.a), b_(factors.b), scratch_size_(factors.scratch),
          a_plan_(a_), b_plan_(b_), input_map_(n_), output_map_(n_) {
        const auto inv_b_a = detail::modular_inverse(b_ % a_, a_);
        const auto inv_a_b = detail::modular_inverse(a_ % b_, b_);
        for (std::size_t n1 = 0; n1 < a_; ++n1) {
            for (std::size_t n2 = 0; n2 < b_; ++n2) {
                const auto i = n1 * b_ + n2;
                const auto first = detail::mul_mod(detail::mul_mod(n1, b_, n_), inv_b_a, n_);
                const auto second = detail::mul_mod(detail::mul_mod(n2, a_, n_), inv_a_b, n_);
                input_map_[i] = detail::add_mod(first, second, n_);
            }
        }
        for (std::size_t k1 = 0; k1 < a_; ++k1) {
            for (std::size_t k2 = 0; k2 < b_; ++k2) {
                output_map_[k1 * b_ + k2] = detail::add_mod(
                    detail::mul_mod(k1, b_, n_), detail::mul_mod(k2, a_, n_), n_);
            }
        }
    }

    void execute(std::span<Complex<T>> data, std::span<Complex<T>> scratch,
                 Direction direction) const {
        if (data.size() != n_ || scratch.size() < scratch_size_)
            throw std::invalid_argument("GoodThomasPlan buffer size mismatch");
        auto matrix = scratch.first(n_);
        const auto max_factor = std::max(a_, b_);
        auto work = scratch.subspan(n_, max_factor);
        auto work_scratch = scratch.subspan(n_ + max_factor, max_factor);
        for (std::size_t i = 0; i < n_; ++i) matrix[i] = data[input_map_[i]];
        for (std::size_t n1 = 0; n1 < a_; ++n1) {
            auto row = matrix.subspan(n1 * b_, b_);
            if (direction == Direction::forward) b_plan_.forward_inplace(row, work_scratch.first(b_));
            else b_plan_.inverse_inplace(row, work_scratch.first(b_));
        }
        for (std::size_t k2 = 0; k2 < b_; ++k2) {
            for (std::size_t n1 = 0; n1 < a_; ++n1) work[n1] = matrix[n1 * b_ + k2];
            if (direction == Direction::forward) a_plan_.forward_inplace(work.first(a_), work_scratch.first(a_));
            else a_plan_.inverse_inplace(work.first(a_), work_scratch.first(a_));
            for (std::size_t k1 = 0; k1 < a_; ++k1) matrix[k1 * b_ + k2] = work[k1];
        }
        for (std::size_t i = 0; i < n_; ++i) data[output_map_[i]] = matrix[i];
    }
    void copy_checked(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                      std::span<Complex<T>> scratch) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < scratch_size_)
            throw std::invalid_argument("GoodThomasPlan buffer size mismatch");
        std::copy(input.begin(), input.end(), output.begin());
    }

    std::size_t n_{}, a_{}, b_{}, scratch_size_{};
    MixedRadixPlan<T> a_plan_, b_plan_;
    std::vector<std::size_t> input_map_, output_map_;
};
GoodThomasPlan(std::size_t) -> GoodThomasPlan<double>;
GoodThomasPlan(std::size_t, std::size_t, std::size_t) -> GoodThomasPlan<double>;

template <Scalar T = double>
class BluesteinPlan {
public:
    explicit BluesteinPlan(std::size_t n)
        : n_(n), m_(workspace(n)), convolution_plan_(m_), chirp_(n), kernel_spectrum_(m_) {
        for (std::size_t k = 0; k < n_; ++k) {
            const auto c = forward_chirp(k, n_);
            chirp_[k] = c;
            const auto b = std::conj(c);
            kernel_spectrum_[k] = b;
            if (k != 0) kernel_spectrum_[m_ - k] = b;
        }
        convolution_plan_.forward_inplace(kernel_spectrum_);
    }
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t convolution_size() const noexcept { return m_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return m_; }
    void forward(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const { execute(input, output, scratch, Direction::forward); }
    void inverse(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const { execute(input, output, scratch, Direction::inverse); }
private:
    [[nodiscard]] static std::size_t workspace(std::size_t n) {
        if (n == 0) throw std::invalid_argument("BluesteinPlan requires length >= 1");
        if (n == 1) return 1;
        if (n > (std::numeric_limits<std::size_t>::max() / 2) + 1)
            throw std::length_error("BluesteinPlan workspace overflow");
        return detail::next_power_of_two(n + (n - 1));
    }
    [[nodiscard]] static Complex<T> forward_chirp(std::size_t k, std::size_t n) {
        const T phase = detail::chirp_phase_ratio<T>(k, n);
        return detail::root<T>(-std::numbers::pi_v<T> * phase);
    }
    void execute(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch, Direction direction) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < m_)
            throw std::invalid_argument("BluesteinPlan buffer size mismatch");
        auto work = scratch.first(m_);
        std::fill(work.begin(), work.end(), Complex<T>{});
        for (std::size_t k = 0; k < n_; ++k) {
            const auto x = direction == Direction::inverse ? std::conj(input[k]) : input[k];
            work[k] = x * chirp_[k];
        }
        convolution_plan_.forward_inplace(work);
        for (std::size_t k = 0; k < m_; ++k) work[k] *= kernel_spectrum_[k];
        convolution_plan_.inverse_inplace(work);
        if (direction == Direction::forward) {
            for (std::size_t k = 0; k < n_; ++k) output[k] = work[k] * chirp_[k];
        } else {
            const T scale = T{1} / static_cast<T>(n_);
            for (std::size_t k = 0; k < n_; ++k) output[k] = std::conj(work[k] * chirp_[k]) * scale;
        }
    }
    std::size_t n_{}, m_{};
    Radix2Plan<T> convolution_plan_;
    Vector<T> chirp_, kernel_spectrum_;
};
BluesteinPlan(std::size_t) -> BluesteinPlan<double>;

template <Scalar T = double>
class RaderPlan {
public:
    explicit RaderPlan(std::size_t n)
        : n_(n), l_(n > 0 ? n - 1 : 0), m_(workspace_checked(n)), direct_cyclic_(m_ == l_),
          convolution_plan_(m_), output_permutation_(l_), input_permutation_(l_), kernel_spectrum_(m_) {
        const auto g = detail::primitive_root_prime(n_);
        output_permutation_[0] = 1;
        for (std::size_t q = 1; q < l_; ++q)
            output_permutation_[q] = detail::mul_mod(output_permutation_[q - 1], g, n_);
        for (std::size_t q = 0; q < l_; ++q) input_permutation_[q] = output_permutation_[(l_ - q) % l_];
        for (std::size_t q = 0; q < l_; ++q)
            kernel_spectrum_[q] = detail::root<T>(
                -T{2} * std::numbers::pi_v<T> * static_cast<T>(output_permutation_[q]) / static_cast<T>(n_));
        convolution_plan_.forward_inplace(kernel_spectrum_);
    }
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] std::size_t cyclic_size() const noexcept { return l_; }
    [[nodiscard]] std::size_t convolution_size() const noexcept { return m_; }
    [[nodiscard]] bool direct_cyclic_fft() const noexcept { return direct_cyclic_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return m_; }
    void forward(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const { execute(input, output, scratch, Direction::forward); }
    void inverse(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch) const { execute(input, output, scratch, Direction::inverse); }
private:
    [[nodiscard]] static std::size_t workspace_checked(std::size_t n) {
        if (n < 3 || !detail::is_prime(n)) throw std::invalid_argument("RaderPlan requires prime length >= 3");
        const auto l = n - 1;
        if (detail::is_power_of_two(l)) return l;
        if (l > (std::numeric_limits<std::size_t>::max() / 2) + 1)
            throw std::length_error("RaderPlan workspace overflow");
        return detail::next_power_of_two(l + (l - 1));
    }
    void execute(std::span<const Complex<T>> input, std::span<Complex<T>> output,
                 std::span<Complex<T>> scratch, Direction direction) const {
        if (input.size() != n_ || output.size() != n_ || scratch.size() < m_)
            throw std::invalid_argument("RaderPlan buffer size mismatch");
        auto work = scratch.first(m_);
        std::fill(work.begin(), work.end(), Complex<T>{});
        Complex<T> dc{};
        for (std::size_t i = 0; i < n_; ++i) dc += direction == Direction::inverse ? std::conj(input[i]) : input[i];
        for (std::size_t q = 0; q < l_; ++q)
            work[q] = direction == Direction::inverse ? std::conj(input[input_permutation_[q]]) : input[input_permutation_[q]];
        convolution_plan_.forward_inplace(work);
        for (std::size_t i = 0; i < m_; ++i) work[i] *= kernel_spectrum_[i];
        convolution_plan_.inverse_inplace(work);
        const auto x0 = direction == Direction::inverse ? std::conj(input[0]) : input[0];
        if (direction == Direction::forward) {
            output[0] = dc;
            for (std::size_t q = 0; q < l_; ++q) {
                auto c = work[q];
                if (!direct_cyclic_ && q < l_ - 1) c += work[q + l_];
                output[output_permutation_[q]] = x0 + c;
            }
        } else {
            const T scale = T{1} / static_cast<T>(n_);
            output[0] = std::conj(dc) * scale;
            for (std::size_t q = 0; q < l_; ++q) {
                auto c = work[q];
                if (!direct_cyclic_ && q < l_ - 1) c += work[q + l_];
                output[output_permutation_[q]] = std::conj(x0 + c) * scale;
            }
        }
    }
    std::size_t n_{}, l_{}, m_{};
    bool direct_cyclic_{};
    Radix2Plan<T> convolution_plan_;
    std::vector<std::size_t> output_permutation_, input_permutation_;
    Vector<T> kernel_spectrum_;
};
RaderPlan(std::size_t) -> RaderPlan<double>;

}  // namespace signal_processing::fft
