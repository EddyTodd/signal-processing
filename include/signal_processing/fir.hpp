#pragma once

#include "signal_processing/convolution.hpp"
#include "signal_processing/detail/sample.hpp"
#include "signal_processing/fft.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace signal_processing::fir {

// Finite impulse response filter.
// https://en.wikipedia.org/wiki/Finite_impulse_response

namespace detail {

template <signal_processing::detail::Sample T>
inline void require_coefficients(std::span<const T> coefficients) {
    if (coefficients.empty())
        throw std::invalid_argument("FIR filtering requires at least one coefficient");
}

}  // namespace detail

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> direct(std::span<const T> input,
                                            std::span<const T> coefficients) {
    detail::require_coefficients<T>(coefficients);
    std::vector<T> output(input.size(), T{});
    for (std::size_t n = 0; n < input.size(); ++n) {
        T sum{};
        const std::size_t taps = std::min(coefficients.size(), n + 1);
        for (std::size_t k = 0; k < taps; ++k) sum += coefficients[k] * input[n - k];
        output[n] = sum;
    }
    return output;
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> full_direct(std::span<const T> input,
                                                 std::span<const T> coefficients) {
    detail::require_coefficients<T>(coefficients);
    return convolution::direct<T>(input, coefficients);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> fft(std::span<const T> input,
                                         std::span<const T> coefficients) {
    detail::require_coefficients<T>(coefficients);
    if (input.empty()) return {};
    auto output = convolution::fft<T>(input, coefficients);
    output.resize(input.size());
    return output;
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> overlap_save(std::span<const T> input,
                                                  std::span<const T> coefficients,
                                                  std::size_t transform_size) {
    detail::require_coefficients<T>(coefficients);
    if (input.empty()) return {};
    auto output = convolution::overlap_save<T>(input, coefficients, transform_size);
    output.resize(input.size());
    return output;
}

template <fft::Scalar T>
[[nodiscard]] inline bool symmetric_coefficients(std::span<const T> coefficients,
                                                  T tolerance = T{8} *
                                                      std::numeric_limits<T>::epsilon()) {
    if (!(tolerance >= T{0}) || !std::isfinite(tolerance))
        throw std::invalid_argument("symmetry tolerance must be finite and nonnegative");
    if (coefficients.empty()) return false;
    for (std::size_t k = 0; k < coefficients.size() / 2; ++k) {
        const T a = coefficients[k];
        const T b = coefficients[coefficients.size() - 1 - k];
        const T scale = T{1} + std::max(std::abs(a), std::abs(b));
        if (std::abs(a - b) > tolerance * scale) return false;
    }
    return true;
}

template <fft::Scalar T>
[[nodiscard]] inline bool antisymmetric_coefficients(std::span<const T> coefficients,
                                                      T tolerance = T{8} *
                                                          std::numeric_limits<T>::epsilon()) {
    if (!(tolerance >= T{0}) || !std::isfinite(tolerance))
        throw std::invalid_argument("antisymmetry tolerance must be finite and nonnegative");
    if (coefficients.empty()) return false;
    for (std::size_t k = 0; k < coefficients.size() / 2; ++k) {
        const T a = coefficients[k];
        const T b = coefficients[coefficients.size() - 1 - k];
        const T scale = T{1} + std::max(std::abs(a), std::abs(b));
        if (std::abs(a + b) > tolerance * scale) return false;
    }
    if ((coefficients.size() & 1U) != 0U &&
        std::abs(coefficients[coefficients.size() / 2]) > tolerance) return false;
    return true;
}

// Linear-phase FIR evaluation exploiting h[k] = h[M-1-k].
// https://en.wikipedia.org/wiki/Finite_impulse_response#Linear_phase

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> symmetric_direct(std::span<const T> input,
                                                      std::span<const T> coefficients) {
    if (!symmetric_coefficients<T>(coefficients))
        throw std::invalid_argument("symmetric_direct requires nonempty symmetric coefficients");
    std::vector<T> output(input.size(), T{});
    const std::size_t m = coefficients.size();
    const std::size_t pairs = m / 2;
    for (std::size_t n = 0; n < input.size(); ++n) {
        T sum{};
        for (std::size_t k = 0; k < pairs; ++k) {
            const std::size_t mirror = m - 1 - k;
            const T first = n >= k ? input[n - k] : T{};
            const T second = n >= mirror ? input[n - mirror] : T{};
            sum += coefficients[k] * (first + second);
        }
        if ((m & 1U) != 0U && n >= pairs) sum += coefficients[pairs] * input[n - pairs];
        output[n] = sum;
    }
    return output;
}

// Linear-phase FIR evaluation exploiting h[k] = -h[M-1-k].
template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> antisymmetric_direct(std::span<const T> input,
                                                          std::span<const T> coefficients) {
    if (!antisymmetric_coefficients<T>(coefficients))
        throw std::invalid_argument("antisymmetric_direct requires nonempty antisymmetric coefficients");
    std::vector<T> output(input.size(), T{});
    const std::size_t m = coefficients.size();
    const std::size_t pairs = m / 2;
    for (std::size_t n = 0; n < input.size(); ++n) {
        T sum{};
        for (std::size_t k = 0; k < pairs; ++k) {
            const std::size_t mirror = m - 1 - k;
            const T first = n >= k ? input[n - k] : T{};
            const T second = n >= mirror ? input[n - mirror] : T{};
            sum += coefficients[k] * (first - second);
        }
        output[n] = sum;
    }
    return output;
}

template <signal_processing::detail::Sample T>
class Filter {
public:
    explicit Filter(std::vector<T> coefficients)
        : coefficients_(std::move(coefficients)), delay_(coefficients_.size(), T{}) {
        if (coefficients_.empty())
            throw std::invalid_argument("FIR filter requires at least one coefficient");
    }

    [[nodiscard]] T process(T sample) noexcept {
        delay_[write_] = sample;
        T output{};
        std::size_t index = write_;
        for (const T& coefficient : coefficients_) {
            output += coefficient * delay_[index];
            index = index == 0 ? delay_.size() - 1 : index - 1;
        }
        write_ = write_ + 1 == delay_.size() ? 0 : write_ + 1;
        return output;
    }

    void process(std::span<const T> input, std::span<T> output) {
        if (output.size() < input.size())
            throw std::invalid_argument("FIR output span is too small");
        for (std::size_t i = 0; i < input.size(); ++i) output[i] = process(input[i]);
    }

    [[nodiscard]] std::vector<T> process(std::span<const T> input) {
        std::vector<T> output(input.size());
        process(input, output);
        return output;
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

template <signal_processing::detail::Sample T>
class OverlapSaveFilter {
public:
    OverlapSaveFilter(std::span<const T> coefficients, std::size_t transform_size)
        : plan_(coefficients, transform_size) {}

    [[nodiscard]] std::size_t block_size() const noexcept { return plan_.block_size(); }
    [[nodiscard]] std::size_t transform_size() const noexcept { return plan_.transform_size(); }
    [[nodiscard]] std::size_t coefficient_count() const noexcept { return plan_.kernel_size(); }

    void reset() { plan_.reset(); }

    void process_block(std::span<const T> input, std::span<T> output) {
        plan_.process_block(input, output);
    }

private:
    convolution::OverlapSavePlan<T> plan_;
};

}  // namespace signal_processing::fir
