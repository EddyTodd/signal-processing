#pragma once

#include <algorithm>
#include <array>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace signal_processing::iir {

// Infinite impulse response filter.
// https://en.wikipedia.org/wiki/Infinite_impulse_response

template <std::floating_point T>
class DirectFormI {
public:
    DirectFormI(std::vector<T> numerator, std::vector<T> denominator)
        : b_(std::move(numerator)), a_(std::move(denominator)),
          x_(b_.size(), T{}), y_(a_.size() > 1 ? a_.size() - 1 : 0, T{}) {
        if (b_.empty() || a_.empty() || a_.front() == T{}) {
            throw std::invalid_argument("IIR coefficients require non-empty b and non-zero a[0]");
        }
        const T a0 = a_.front();
        if (a0 != T{1}) {
            for (auto& coefficient : b_) {
                coefficient /= a0;
            }
            for (auto& coefficient : a_) {
                coefficient /= a0;
            }
        }
    }

    [[nodiscard]] T process(T sample) noexcept {
        for (std::size_t i = x_.size(); i-- > 1;) {
            x_[i] = x_[i - 1];
        }
        x_[0] = sample;

        T output{};
        for (std::size_t i = 0; i < b_.size(); ++i) {
            output += b_[i] * x_[i];
        }
        for (std::size_t i = 1; i < a_.size(); ++i) {
            output -= a_[i] * y_[i - 1];
        }

        for (std::size_t i = y_.size(); i-- > 1;) {
            y_[i] = y_[i - 1];
        }
        if (!y_.empty()) {
            y_[0] = output;
        }
        return output;
    }

    void reset() noexcept {
        std::fill(x_.begin(), x_.end(), T{});
        std::fill(y_.begin(), y_.end(), T{});
    }

private:
    std::vector<T> b_;
    std::vector<T> a_;
    std::vector<T> x_;
    std::vector<T> y_;
};

template <std::floating_point T>
class BiquadTransposedDirectFormII {
public:
    BiquadTransposedDirectFormII(T b0, T b1, T b2, T a0, T a1, T a2) {
        if (a0 == T{}) {
            throw std::invalid_argument("biquad requires non-zero a0");
        }
        b0_ = b0 / a0;
        b1_ = b1 / a0;
        b2_ = b2 / a0;
        a1_ = a1 / a0;
        a2_ = a2 / a0;
    }

    [[nodiscard]] T process(T sample) noexcept {
        const T output = b0_ * sample + s1_;
        s1_ = b1_ * sample - a1_ * output + s2_;
        s2_ = b2_ * sample - a2_ * output;
        return output;
    }

    void reset() noexcept {
        s1_ = T{};
        s2_ = T{};
    }

private:
    T b0_{};
    T b1_{};
    T b2_{};
    T a1_{};
    T a2_{};
    T s1_{};
    T s2_{};
};

}  // namespace signal_processing::iir
