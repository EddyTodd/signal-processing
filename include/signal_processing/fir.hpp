#pragma once

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace signal_processing::fir {

// Finite impulse response filter.
// https://en.wikipedia.org/wiki/Finite_impulse_response

template <std::floating_point T>
[[nodiscard]] inline std::vector<T> direct(std::span<const T> input,
                                            std::span<const T> coefficients) {
    std::vector<T> output(input.size(), T{});
    for (std::size_t n = 0; n < input.size(); ++n) {
        T sum{};
        const std::size_t taps = std::min(coefficients.size(), n + 1);
        for (std::size_t k = 0; k < taps; ++k) {
            sum += coefficients[k] * input[n - k];
        }
        output[n] = sum;
    }
    return output;
}

template <std::floating_point T>
class Filter {
public:
    explicit Filter(std::vector<T> coefficients)
        : coefficients_(std::move(coefficients)), delay_(coefficients_.size(), T{}) {
        if (coefficients_.empty()) {
            throw std::invalid_argument("FIR filter requires at least one coefficient");
        }
    }

    [[nodiscard]] T process(T sample) noexcept {
        delay_[write_] = sample;
        T output{};
        std::size_t index = write_;
        for (const T coefficient : coefficients_) {
            output += coefficient * delay_[index];
            index = index == 0 ? delay_.size() - 1 : index - 1;
        }
        write_ = (write_ + 1) % delay_.size();
        return output;
    }

    void process(std::span<const T> input, std::span<T> output) {
        if (output.size() < input.size()) {
            throw std::invalid_argument("FIR output span is too small");
        }
        for (std::size_t i = 0; i < input.size(); ++i) {
            output[i] = process(input[i]);
        }
    }

    void reset() noexcept {
        std::fill(delay_.begin(), delay_.end(), T{});
        write_ = 0;
    }

    [[nodiscard]] std::span<const T> coefficients() const noexcept { return coefficients_; }

private:
    std::vector<T> coefficients_;
    std::vector<T> delay_;
    std::size_t write_{0};
};

}  // namespace signal_processing::fir
