#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace signal_processing::iir {

// Infinite impulse response filter.
// https://en.wikipedia.org/wiki/Infinite_impulse_response

template <typename T>
concept Scalar = std::same_as<T, float> || std::same_as<T, double>;

template <Scalar T>
struct Coefficients {
    std::vector<T> numerator;
    std::vector<T> denominator;
};

namespace detail {

template <Scalar T>
void require_finite(std::span<const T> values, const char* message) {
    if (!std::all_of(values.begin(), values.end(), [](T value) { return std::isfinite(value); }))
        throw std::invalid_argument(message);
}

template <Scalar T>
void require_finite(T value, const char* message) {
    if (!std::isfinite(value)) throw std::invalid_argument(message);
}

}  // namespace detail

template <Scalar T>
[[nodiscard]] inline Coefficients<T> normalize(std::vector<T> numerator,
                                                std::vector<T> denominator) {
    if (numerator.empty() || denominator.empty())
        throw std::invalid_argument("IIR coefficients require nonempty numerator and denominator");
    detail::require_finite<T>(numerator, "IIR numerator coefficients must be finite");
    detail::require_finite<T>(denominator, "IIR denominator coefficients must be finite");
    if (denominator.front() == T{})
        throw std::invalid_argument("IIR denominator requires nonzero a[0]");
    const T a0 = denominator.front();
    for (auto& value : numerator) value /= a0;
    for (auto& value : denominator) value /= a0;
    return {std::move(numerator), std::move(denominator)};
}

template <Scalar T>
class DirectFormI {
public:
    DirectFormI(std::vector<T> numerator, std::vector<T> denominator) {
        auto coefficients = normalize<T>(std::move(numerator), std::move(denominator));
        b_ = std::move(coefficients.numerator);
        a_ = std::move(coefficients.denominator);
        x_.assign(b_.size(), T{});
        y_.assign(a_.size() > 1 ? a_.size() - 1 : 0, T{});
    }

    [[nodiscard]] T process(T sample) noexcept {
        for (std::size_t i = x_.size(); i-- > 1;) x_[i] = x_[i - 1];
        x_[0] = sample;

        T output{};
        for (std::size_t i = 0; i < b_.size(); ++i) output += b_[i] * x_[i];
        for (std::size_t i = 1; i < a_.size(); ++i) output -= a_[i] * y_[i - 1];

        for (std::size_t i = y_.size(); i-- > 1;) y_[i] = y_[i - 1];
        if (!y_.empty()) y_[0] = output;
        return output;
    }

    void process(std::span<const T> input, std::span<T> output) {
        if (output.size() < input.size())
            throw std::invalid_argument("IIR output span is too small");
        for (std::size_t i = 0; i < input.size(); ++i) output[i] = process(input[i]);
    }

    void reset() noexcept {
        std::fill(x_.begin(), x_.end(), T{});
        std::fill(y_.begin(), y_.end(), T{});
    }

private:
    std::vector<T> b_, a_, x_, y_;
};

template <Scalar T>
class DirectFormII {
public:
    DirectFormII(std::vector<T> numerator, std::vector<T> denominator) {
        auto coefficients = normalize<T>(std::move(numerator), std::move(denominator));
        const std::size_t order =
            std::max(coefficients.numerator.size(), coefficients.denominator.size()) - 1;
        b_.assign(order + 1, T{});
        a_.assign(order + 1, T{});
        std::copy(coefficients.numerator.begin(), coefficients.numerator.end(), b_.begin());
        std::copy(coefficients.denominator.begin(), coefficients.denominator.end(), a_.begin());
        state_.assign(order, T{});
    }

    [[nodiscard]] T process(T sample) noexcept {
        T w = sample;
        for (std::size_t i = 1; i < a_.size(); ++i) w -= a_[i] * state_[i - 1];
        T output = b_[0] * w;
        for (std::size_t i = 1; i < b_.size(); ++i) output += b_[i] * state_[i - 1];
        for (std::size_t i = state_.size(); i-- > 1;) state_[i] = state_[i - 1];
        if (!state_.empty()) state_[0] = w;
        return output;
    }

    void process(std::span<const T> input, std::span<T> output) {
        if (output.size() < input.size())
            throw std::invalid_argument("IIR output span is too small");
        for (std::size_t i = 0; i < input.size(); ++i) output[i] = process(input[i]);
    }

    void reset() noexcept { std::fill(state_.begin(), state_.end(), T{}); }

private:
    std::vector<T> b_, a_, state_;
};

// Transpose of the two-delay-line Direct Form I signal-flow graph.
template <Scalar T>
class TransposedDirectFormI {
public:
    TransposedDirectFormI(std::vector<T> numerator, std::vector<T> denominator) {
        auto coefficients = normalize<T>(std::move(numerator), std::move(denominator));
        b_ = std::move(coefficients.numerator);
        a_ = std::move(coefficients.denominator);
        feedforward_.assign(b_.size() > 1 ? b_.size() - 1 : 0, T{});
        feedback_.assign(a_.size() > 1 ? a_.size() - 1 : 0, T{});
    }

    [[nodiscard]] T process(T sample) noexcept {
        const T internal = sample + (feedback_.empty() ? T{} : feedback_[0]);
        const T output = b_[0] * internal + (feedforward_.empty() ? T{} : feedforward_[0]);

        for (std::size_t i = 0; i < feedforward_.size(); ++i) {
            const T next = i + 1 < feedforward_.size() ? feedforward_[i + 1] : T{};
            feedforward_[i] = next + b_[i + 1] * internal;
        }
        for (std::size_t i = 0; i < feedback_.size(); ++i) {
            const T next = i + 1 < feedback_.size() ? feedback_[i + 1] : T{};
            feedback_[i] = next - a_[i + 1] * internal;
        }
        return output;
    }

    void process(std::span<const T> input, std::span<T> output) {
        if (output.size() < input.size())
            throw std::invalid_argument("IIR output span is too small");
        for (std::size_t i = 0; i < input.size(); ++i) output[i] = process(input[i]);
    }

    void reset() noexcept {
        std::fill(feedforward_.begin(), feedforward_.end(), T{});
        std::fill(feedback_.begin(), feedback_.end(), T{});
    }

private:
    std::vector<T> b_, a_, feedforward_, feedback_;
};

template <Scalar T>
class TransposedDirectFormII {
public:
    TransposedDirectFormII(std::vector<T> numerator, std::vector<T> denominator) {
        auto coefficients = normalize<T>(std::move(numerator), std::move(denominator));
        const std::size_t order =
            std::max(coefficients.numerator.size(), coefficients.denominator.size()) - 1;
        b_.assign(order + 1, T{});
        a_.assign(order + 1, T{});
        std::copy(coefficients.numerator.begin(), coefficients.numerator.end(), b_.begin());
        std::copy(coefficients.denominator.begin(), coefficients.denominator.end(), a_.begin());
        state_.assign(order, T{});
    }

    [[nodiscard]] T process(T sample) noexcept {
        if (state_.empty()) return b_[0] * sample;
        const T output = b_[0] * sample + state_[0];
        for (std::size_t i = 0; i + 1 < state_.size(); ++i)
            state_[i] = state_[i + 1] + b_[i + 1] * sample - a_[i + 1] * output;
        const std::size_t last = state_.size() - 1;
        state_[last] = b_[last + 1] * sample - a_[last + 1] * output;
        return output;
    }

    void process(std::span<const T> input, std::span<T> output) {
        if (output.size() < input.size())
            throw std::invalid_argument("IIR output span is too small");
        for (std::size_t i = 0; i < input.size(); ++i) output[i] = process(input[i]);
    }

    void reset() noexcept { std::fill(state_.begin(), state_.end(), T{}); }

private:
    std::vector<T> b_, a_, state_;
};

template <Scalar T>
struct FirstOrderCoefficients {
    T b0{};
    T b1{};
    T a0{1};
    T a1{};
};

template <Scalar T>
class FirstOrderTransposedDirectFormII {
public:
    explicit FirstOrderTransposedDirectFormII(FirstOrderCoefficients<T> c) {
        if (!std::isfinite(c.b0) || !std::isfinite(c.b1) || !std::isfinite(c.a0) ||
            !std::isfinite(c.a1))
            throw std::invalid_argument("first-order section coefficients must be finite");
        if (c.a0 == T{}) throw std::invalid_argument("first-order section requires nonzero a0");
        b0_ = c.b0 / c.a0;
        b1_ = c.b1 / c.a0;
        a1_ = c.a1 / c.a0;
    }

    [[nodiscard]] T process(T sample) noexcept {
        const T output = b0_ * sample + state_;
        state_ = b1_ * sample - a1_ * output;
        return output;
    }

    void reset() noexcept { state_ = T{}; }

private:
    T b0_{}, b1_{}, a1_{}, state_{};
};

// Digital biquad filter.
// https://en.wikipedia.org/wiki/Digital_biquad_filter
template <Scalar T>
struct BiquadCoefficients {
    T b0{};
    T b1{};
    T b2{};
    T a0{1};
    T a1{};
    T a2{};
};

namespace detail {

template <Scalar T>
void validate_biquad(BiquadCoefficients<T> c) {
    if (!std::isfinite(c.b0) || !std::isfinite(c.b1) || !std::isfinite(c.b2) ||
        !std::isfinite(c.a0) || !std::isfinite(c.a1) || !std::isfinite(c.a2))
        throw std::invalid_argument("biquad coefficients must be finite");
    if (c.a0 == T{}) throw std::invalid_argument("biquad requires nonzero a0");
}

}  // namespace detail

template <Scalar T>
class BiquadDirectFormI {
public:
    explicit BiquadDirectFormI(BiquadCoefficients<T> c) {
        detail::validate_biquad(c);
        b0_ = c.b0 / c.a0;
        b1_ = c.b1 / c.a0;
        b2_ = c.b2 / c.a0;
        a1_ = c.a1 / c.a0;
        a2_ = c.a2 / c.a0;
    }

    [[nodiscard]] T process(T x) noexcept {
        const T y = b0_ * x + b1_ * x1_ + b2_ * x2_ - a1_ * y1_ - a2_ * y2_;
        x2_ = x1_;
        x1_ = x;
        y2_ = y1_;
        y1_ = y;
        return y;
    }

    void reset() noexcept { x1_ = x2_ = y1_ = y2_ = T{}; }

private:
    T b0_{}, b1_{}, b2_{}, a1_{}, a2_{}, x1_{}, x2_{}, y1_{}, y2_{};
};

template <Scalar T>
class BiquadDirectFormII {
public:
    explicit BiquadDirectFormII(BiquadCoefficients<T> c) {
        detail::validate_biquad(c);
        b0_ = c.b0 / c.a0;
        b1_ = c.b1 / c.a0;
        b2_ = c.b2 / c.a0;
        a1_ = c.a1 / c.a0;
        a2_ = c.a2 / c.a0;
    }

    [[nodiscard]] T process(T x) noexcept {
        const T w = x - a1_ * w1_ - a2_ * w2_;
        const T y = b0_ * w + b1_ * w1_ + b2_ * w2_;
        w2_ = w1_;
        w1_ = w;
        return y;
    }

    void reset() noexcept { w1_ = w2_ = T{}; }

private:
    T b0_{}, b1_{}, b2_{}, a1_{}, a2_{}, w1_{}, w2_{};
};

template <Scalar T>
class BiquadTransposedDirectFormII {
public:
    BiquadTransposedDirectFormII(T b0, T b1, T b2, T a0, T a1, T a2)
        : BiquadTransposedDirectFormII(BiquadCoefficients<T>{b0, b1, b2, a0, a1, a2}) {}

    explicit BiquadTransposedDirectFormII(BiquadCoefficients<T> c) {
        detail::validate_biquad(c);
        b0_ = c.b0 / c.a0;
        b1_ = c.b1 / c.a0;
        b2_ = c.b2 / c.a0;
        a1_ = c.a1 / c.a0;
        a2_ = c.a2 / c.a0;
    }

    [[nodiscard]] T process(T sample) noexcept {
        const T output = b0_ * sample + s1_;
        s1_ = b1_ * sample - a1_ * output + s2_;
        s2_ = b2_ * sample - a2_ * output;
        return output;
    }

    void reset() noexcept { s1_ = s2_ = T{}; }

private:
    T b0_{}, b1_{}, b2_{}, a1_{}, a2_{}, s1_{}, s2_{};
};

template <Scalar T>
class SosCascade {
public:
    explicit SosCascade(std::vector<BiquadCoefficients<T>> sections, T gain = T{1})
        : coefficients_(std::move(sections)), gain_(gain) {
        detail::require_finite<T>(gain_, "SOS gain must be finite");
        filters_.reserve(coefficients_.size());
        for (const auto& c : coefficients_) filters_.emplace_back(c);
    }

    [[nodiscard]] T process(T sample) noexcept {
        T value = gain_ * sample;
        for (auto& section : filters_) value = section.process(value);
        return value;
    }

    void process(std::span<const T> input, std::span<T> output) {
        if (output.size() < input.size())
            throw std::invalid_argument("SOS output span is too small");
        for (std::size_t i = 0; i < input.size(); ++i) output[i] = process(input[i]);
    }

    void reset() noexcept {
        for (auto& section : filters_) section.reset();
    }

    [[nodiscard]] std::span<const BiquadCoefficients<T>> sections() const noexcept {
        return coefficients_;
    }
    [[nodiscard]] T gain() const noexcept { return gain_; }

private:
    std::vector<BiquadCoefficients<T>> coefficients_;
    T gain_{1};
    std::vector<BiquadTransposedDirectFormII<T>> filters_;
};

}  // namespace signal_processing::iir
