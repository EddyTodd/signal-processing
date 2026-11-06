#pragma once

#include "signal_processing/convolution.hpp"
#include "signal_processing/detail/sample.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <span>
#include <vector>

namespace signal_processing::correlation {

namespace detail {

// Range energies must not use prefix[right] - prefix[left]:
// cancellation can erase a weak signal at the end of a loud stream.
// A pairwise range-sum tree evaluates only samples inside the overlap,
// with O(N) storage and O(log N) per correlation lag.
template <signal_processing::detail::Sample T>
class EnergyTree {
public:
    explicit EnergyTree(std::span<const T> input)
        : count_(input.size()), sums_(checked_capacity(count_), 0.L) {
        for (std::size_t i = 0; i < count_; ++i) {
            if constexpr (signal_processing::detail::SampleTraits<T>::is_complex) {
                const auto real = static_cast<long double>(input[i].real());
                const auto imag = static_cast<long double>(input[i].imag());
                sums_[count_ + i] = real * real + imag * imag;
            } else {
                const auto value = static_cast<long double>(input[i]);
                sums_[count_ + i] = value * value;
            }
        }
        for (std::size_t i = count_; i > 1; --i)
            sums_[i - 1] = sums_[2 * (i - 1)] + sums_[2 * (i - 1) + 1];
    }

    [[nodiscard]] long double sum(std::size_t first, std::size_t last) const {
        // The caller computes mathematically valid overlapping intervals.
        long double result{};
        first += count_;
        last += count_;
        while (first < last) {
            if ((first & 1U) != 0U) result += sums_[first++];
            if ((last & 1U) != 0U) result += sums_[--last];
            first >>= 1U;
            last >>= 1U;
        }
        return result;
    }

private:
    [[nodiscard]] static std::size_t checked_capacity(std::size_t count) {
        if (count > std::numeric_limits<std::size_t>::max() / 2)
            throw std::length_error("normalized correlation energy workspace overflows");
        return 2 * count;
    }

    std::size_t count_{};
    std::vector<long double> sums_;
};

template <signal_processing::detail::Sample T>
[[nodiscard]] std::vector<T> normalize(std::span<const T> lhs, std::span<const T> rhs,
                                       std::span<const T> raw, bool from_fft = false) {
    std::vector<T> output(raw.size(), T{});
    if (raw.empty()) return output;
    const EnergyTree<T> lhs_energy(lhs);
    const EnergyTree<T> rhs_energy(rhs);
    using Scalar = signal_processing::detail::scalar_t<T>;
    // FFT convolution has a nonzero absolute roundoff floor: a tiny
    // overlap can give a large *relative* normalized error even though
    // the unnormalized Fourier result satisfies a normal absolute bound.
    // Directly recompute only ill-conditioned overlaps in long double.
    const long double total_energy_scale =
        std::sqrt(lhs_energy.sum(0, lhs.size())) *
        std::sqrt(rhs_energy.sum(0, rhs.size()));
    const long double fft_roundoff_floor =
        from_fft ? 4096.L *
            static_cast<long double>(std::numeric_limits<Scalar>::epsilon()) *
            (2.L + std::log2(1.L + static_cast<long double>(
                std::max(lhs.size(), rhs.size())))) * total_energy_scale
            : 0.L;
    const auto lhs_size = static_cast<std::ptrdiff_t>(lhs.size());
    const auto rhs_size = static_cast<std::ptrdiff_t>(rhs.size());

    for (std::size_t index = 0; index < raw.size(); ++index) {
        const auto lag = static_cast<std::ptrdiff_t>(index) - (rhs_size - 1);
        const auto i_first = std::max<std::ptrdiff_t>(0, lag);
        const auto i_last = std::min<std::ptrdiff_t>(lhs_size, rhs_size + lag);
        if (i_first >= i_last) continue;
        const auto j_first = i_first - lag;
        const auto j_last = i_last - lag;
        const long double left = lhs_energy.sum(static_cast<std::size_t>(i_first),
                                                static_cast<std::size_t>(i_last));
        const long double right = rhs_energy.sum(static_cast<std::size_t>(j_first),
                                                 static_cast<std::size_t>(j_last));
        // Compute square roots separately to reduce product overflow.
        const long double denominator = std::sqrt(left) * std::sqrt(right);
        if (!(denominator > 0.L)) continue; // Zero-energy overlap maps to zero.
        if (from_fft && denominator <= fft_roundoff_floor) {
            std::complex<long double> numerator{};
            for (auto i = i_first; i < i_last; ++i) {
                const auto j = i - lag;
                const auto x = lhs[static_cast<std::size_t>(i)];
                const auto y = rhs[static_cast<std::size_t>(j)];
                if constexpr (signal_processing::detail::SampleTraits<T>::is_complex) {
                    const std::complex<long double> xl{
                        static_cast<long double>(x.real()),
                        static_cast<long double>(x.imag())};
                    const std::complex<long double> yl{
                        static_cast<long double>(y.real()),
                        static_cast<long double>(y.imag())};
                    numerator += xl * std::conj(yl);
                } else {
                    numerator += std::complex<long double>{
                        static_cast<long double>(x) * static_cast<long double>(y), 0.L};
                }
            }
            if constexpr (signal_processing::detail::SampleTraits<T>::is_complex) {
                output[index] = T{
                    static_cast<Scalar>(numerator.real() / denominator),
                    static_cast<Scalar>(numerator.imag() / denominator)};
            } else {
                output[index] = static_cast<T>(numerator.real() / denominator);
            }
        } else if constexpr (signal_processing::detail::SampleTraits<T>::is_complex) {
            output[index] = T{
                static_cast<Scalar>(static_cast<long double>(raw[index].real()) / denominator),
                static_cast<Scalar>(static_cast<long double>(raw[index].imag()) / denominator)};
        } else {
            output[index] = static_cast<T>(
                static_cast<long double>(raw[index]) / denominator);
        }
    }
    return output;
}

}  // namespace detail

// r_xy[lag] = sum_n x[n] conj(y[n-lag]); output index = lag + rhs.size() - 1.
// https://en.wikipedia.org/wiki/Cross-correlation
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> cross_direct(std::span<const T> lhs,
                                                  std::span<const T> rhs) {
    if (lhs.empty() || rhs.empty()) return {};
    const std::size_t size = convolution::detail::linear_size(lhs.size(), rhs.size());
    std::vector<T> output(size, T{});
    const auto rhs_extent = static_cast<std::ptrdiff_t>(rhs.size());
    for (std::size_t index = 0; index < size; ++index) {
        const auto lag = static_cast<std::ptrdiff_t>(index) - (rhs_extent - 1);
        T sum{};
        for (std::size_t i = 0; i < lhs.size(); ++i) {
            const auto j = static_cast<std::ptrdiff_t>(i) - lag;
            if (j >= 0 && j < rhs_extent)
                sum += lhs[i] * signal_processing::detail::conjugate(
                                    rhs[static_cast<std::size_t>(j)]);
        }
        output[index] = sum;
    }
    return output;
}

// Cross-correlation through convolution with the reversed conjugate of rhs.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> cross_fft(std::span<const T> lhs,
                                               std::span<const T> rhs) {
    if (lhs.empty() || rhs.empty()) return {};
    std::vector<T> reversed(rhs.size());
    for (std::size_t i = 0; i < rhs.size(); ++i)
        reversed[i] = signal_processing::detail::conjugate(rhs[rhs.size() - 1 - i]);
    return convolution::fft<T>(lhs, reversed);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> auto_direct(std::span<const T> input) {
    return cross_direct<T>(input, input);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> auto_fft(std::span<const T> input) {
    return cross_fft<T>(input, input);
}

// Per-lag normalized cross-correlation. Zero-energy overlaps map to zero.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> normalized_direct(std::span<const T> lhs,
                                                       std::span<const T> rhs) {
    const auto raw = cross_direct<T>(lhs, rhs);
    return detail::normalize<T>(lhs, rhs, raw);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> normalized_fft(std::span<const T> lhs,
                                                    std::span<const T> rhs) {
    const auto raw = cross_fft<T>(lhs, rhs);
    return detail::normalize<T>(lhs, rhs, raw, true);
}

// Compatibility names from the initial milestone; both are direct definitions.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> cross(std::span<const T> lhs, std::span<const T> rhs) {
    return cross_direct<T>(lhs, rhs);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> auto_correlation(std::span<const T> input) {
    return auto_direct<T>(input);
}

}  // namespace signal_processing::correlation
