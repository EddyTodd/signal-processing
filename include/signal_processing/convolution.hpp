#pragma once

#include "signal_processing/detail/sample.hpp"
#include "signal_processing/fft.hpp"

#include <algorithm>
#include <complex>
#include <concepts>
#include <cstddef>
#include <limits>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace signal_processing::convolution {

namespace detail {

[[nodiscard]] inline std::size_t linear_size(std::size_t a, std::size_t b) {
    if (a == 0 || b == 0) return 0;
    if (a > std::numeric_limits<std::size_t>::max() - b + 1)
        throw std::length_error("linear convolution size overflows size_t");
    return a + b - 1;
}

[[nodiscard]] inline std::size_t next_power_of_two(std::size_t n) {
    if (n <= 1) return 1;
    const auto high = std::size_t{1} << (std::numeric_limits<std::size_t>::digits - 1);
    if (n > high) throw std::length_error("FFT workspace overflows size_t");
    --n;
    for (std::size_t shift = 1; shift < std::numeric_limits<std::size_t>::digits; shift <<= 1)
        n |= n >> shift;
    return n + 1;
}

template <signal_processing::detail::Sample T>
void fill_complex(std::span<const T> input,
                  std::span<signal_processing::detail::complex_t<T>> output) {
    if (output.size() < input.size()) throw std::invalid_argument("complex workspace too small");
    std::fill(output.begin(), output.end(), signal_processing::detail::complex_t<T>{});
    for (std::size_t i = 0; i < input.size(); ++i)
        output[i] = signal_processing::detail::to_complex(input[i]);
}

template <signal_processing::detail::Sample T>
[[nodiscard]] std::vector<T> extract(
    std::span<const signal_processing::detail::complex_t<T>> input, std::size_t count) {
    if (input.size() < count) throw std::invalid_argument("complex workspace too small");
    std::vector<T> output(count);
    for (std::size_t i = 0; i < count; ++i)
        output[i] = signal_processing::detail::from_complex<T>(input[i]);
    return output;
}

}  // namespace detail

// Linear convolution by the defining O(NM) sum.
// https://en.wikipedia.org/wiki/Convolution
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> direct(std::span<const T> signal,
                                            std::span<const T> kernel) {
    const std::size_t output_size = detail::linear_size(signal.size(), kernel.size());
    if (output_size == 0) return {};
    std::vector<T> output(output_size, T{});
    for (std::size_t i = 0; i < signal.size(); ++i)
        for (std::size_t j = 0; j < kernel.size(); ++j)
            output[i + j] += signal[i] * kernel[j];
    return output;
}

// Equal-length circular convolution by its defining periodic sum.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> circular(std::span<const T> lhs,
                                              std::span<const T> rhs) {
    if (lhs.empty() && rhs.empty()) return {};
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("circular convolution requires equal lengths");
    const std::size_t n = lhs.size();
    std::vector<T> output(n, T{});
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            const std::size_t index = k >= j ? k - j : n - (j - k);
            output[k] += lhs[j] * rhs[index];
        }
    }
    return output;
}

// Zero-padded radix-2 FFT convolution.
// https://en.wikipedia.org/wiki/Convolution_theorem
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> fft(std::span<const T> signal,
                                         std::span<const T> kernel) {
    const std::size_t output_size = detail::linear_size(signal.size(), kernel.size());
    if (output_size == 0) return {};
    const std::size_t n = detail::next_power_of_two(output_size);
    using Scalar = signal_processing::detail::scalar_t<T>;
    using Complex = signal_processing::detail::complex_t<T>;
    std::vector<Complex> a(n), b(n);
    detail::fill_complex<T>(signal, a);
    detail::fill_complex<T>(kernel, b);
    fft::Radix2Plan<Scalar> plan(n);
    plan.forward_inplace(a);
    plan.forward_inplace(b);
    for (std::size_t i = 0; i < n; ++i) a[i] *= b[i];
    plan.inverse_inplace(a);
    return detail::extract<T>(a, output_size);
}

// Arbitrary-length circular convolution using Bluestein's FFT on exactly N samples.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> circular_bluestein(std::span<const T> lhs,
                                                        std::span<const T> rhs) {
    if (lhs.empty() && rhs.empty()) return {};
    if (lhs.size() != rhs.size())
        throw std::invalid_argument("circular convolution requires equal lengths");
    using Scalar = signal_processing::detail::scalar_t<T>;
    using Complex = signal_processing::detail::complex_t<T>;
    std::vector<Complex> a(lhs.size()), b(rhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        a[i] = signal_processing::detail::to_complex(lhs[i]);
        b[i] = signal_processing::detail::to_complex(rhs[i]);
    }
    a = fft::bluestein<Scalar>(a, fft::Direction::forward);
    b = fft::bluestein<Scalar>(b, fft::Direction::forward);
    for (std::size_t i = 0; i < a.size(); ++i) a[i] *= b[i];
    a = fft::bluestein<Scalar>(a, fft::Direction::inverse);
    return detail::extract<T>(a, a.size());
}

// Overlap-add linear convolution. block_size is the number of new input samples per FFT.
// https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> overlap_add(std::span<const T> signal,
                                                 std::span<const T> kernel,
                                                 std::size_t block_size) {
    const std::size_t output_size = detail::linear_size(signal.size(), kernel.size());
    if (output_size == 0) return {};
    if (block_size == 0) throw std::invalid_argument("overlap-add block size must be nonzero");
    if (block_size > std::numeric_limits<std::size_t>::max() - kernel.size() + 1)
        throw std::length_error("overlap-add workspace size overflows size_t");

    const std::size_t transform_size =
        detail::next_power_of_two(block_size + kernel.size() - 1);
    using Scalar = signal_processing::detail::scalar_t<T>;
    using Complex = signal_processing::detail::complex_t<T>;
    fft::Radix2Plan<Scalar> plan(transform_size);
    std::vector<Complex> kernel_spectrum(transform_size), work(transform_size);
    detail::fill_complex<T>(kernel, kernel_spectrum);
    plan.forward_inplace(kernel_spectrum);

    std::vector<T> output(output_size, T{});
    std::size_t offset = 0;
    while (offset < signal.size()) {
        const std::size_t count = std::min(block_size, signal.size() - offset);
        detail::fill_complex<T>(signal.subspan(offset, count), work);
        plan.forward_inplace(work);
        for (std::size_t i = 0; i < transform_size; ++i) work[i] *= kernel_spectrum[i];
        plan.inverse_inplace(work);
        const std::size_t produced = count + kernel.size() - 1;
        for (std::size_t i = 0; i < produced && i < output.size() - offset; ++i)
            output[offset + i] += signal_processing::detail::from_complex<T>(work[i]);
        offset += count;
    }
    return output;
}

template <signal_processing::detail::Sample T>
class OverlapSavePlan {
public:
    OverlapSavePlan(std::span<const T> kernel, std::size_t transform_size)
        : kernel_(kernel.begin(), kernel.end()), transform_size_(transform_size),
          block_size_(checked_block_size(kernel_.size(), transform_size_)),
          plan_(transform_size_), kernel_spectrum_(transform_size_), work_(transform_size_),
          time_block_(transform_size_), overlap_(kernel_.size() - 1, T{}) {
        detail::fill_complex<T>(kernel_, kernel_spectrum_);
        plan_.forward_inplace(kernel_spectrum_);
    }

    [[nodiscard]] std::size_t transform_size() const noexcept { return transform_size_; }
    [[nodiscard]] std::size_t block_size() const noexcept { return block_size_; }
    [[nodiscard]] std::size_t kernel_size() const noexcept { return kernel_.size(); }

    void reset() { std::fill(overlap_.begin(), overlap_.end(), T{}); }

    void process_block(std::span<const T> input, std::span<T> output) {
        if (input.size() != block_size_ || output.size() != block_size_)
            throw std::invalid_argument("OverlapSavePlan block size mismatch");
        const std::size_t overlap_size = overlap_.size();
        std::copy(overlap_.begin(), overlap_.end(), time_block_.begin());
        std::copy(input.begin(), input.end(), std::span<T>(time_block_).subspan(overlap_size).begin());
        if (overlap_size != 0) {
            const auto tail = std::span<const T>(time_block_).last(overlap_size);
            std::copy(tail.begin(), tail.end(), overlap_.begin());
        }
        detail::fill_complex<T>(time_block_, work_);
        plan_.forward_inplace(work_);
        for (std::size_t i = 0; i < transform_size_; ++i) work_[i] *= kernel_spectrum_[i];
        plan_.inverse_inplace(work_);
        for (std::size_t i = 0; i < block_size_; ++i)
            output[i] = signal_processing::detail::from_complex<T>(work_[overlap_size + i]);
    }

private:
    [[nodiscard]] static std::size_t checked_block_size(std::size_t kernel_size,
                                                        std::size_t transform_size) {
        if (kernel_size == 0)
            throw std::invalid_argument("OverlapSavePlan requires a nonempty kernel");
        if (transform_size == 0 || (transform_size & (transform_size - 1)) != 0)
            throw std::invalid_argument("OverlapSavePlan transform size must be a power of two");
        if (transform_size < kernel_size)
            throw std::invalid_argument("OverlapSavePlan transform size must cover the kernel");
        return transform_size - kernel_size + 1;
    }

    using Scalar = signal_processing::detail::scalar_t<T>;
    using Complex = signal_processing::detail::complex_t<T>;
    std::vector<T> kernel_;
    std::size_t transform_size_{};
    std::size_t block_size_{};
    fft::Radix2Plan<Scalar> plan_;
    std::vector<Complex> kernel_spectrum_;
    std::vector<Complex> work_;
    std::vector<T> time_block_;
    std::vector<T> overlap_;
};

// Overlap-save linear convolution. transform_size is the power-of-two FFT length.
// https://en.wikipedia.org/wiki/Overlap%E2%80%93save_method
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> overlap_save(std::span<const T> signal,
                                                  std::span<const T> kernel,
                                                  std::size_t transform_size) {
    const std::size_t output_size = detail::linear_size(signal.size(), kernel.size());
    if (output_size == 0) return {};
    OverlapSavePlan<T> plan(kernel, transform_size);
    const std::size_t block_size = plan.block_size();
    std::vector<T> block(block_size, T{}), produced(block_size, T{}), output;
    output.reserve(output_size);

    std::size_t consumed = 0;
    while (output.size() < output_size) {
        std::fill(block.begin(), block.end(), T{});
        const std::size_t count = std::min(block_size, signal.size() - std::min(consumed, signal.size()));
        if (count != 0) {
            const auto source = signal.subspan(consumed, count);
            std::copy(source.begin(), source.end(), block.begin());
        }
        consumed += count;
        plan.process_block(block, produced);
        const std::size_t keep = std::min(block_size, output_size - output.size());
        const auto valid = std::span<const T>(produced).first(keep);
        output.insert(output.end(), valid.begin(), valid.end());
    }
    return output;
}

template <signal_processing::detail::Sample T>
class StreamingDirect {
public:
    explicit StreamingDirect(std::span<const T> kernel)
        : kernel_(kernel.begin(), kernel.end()), history_(kernel_.size(), T{}) {
        if (kernel_.empty())
            throw std::invalid_argument("StreamingDirect requires a nonempty kernel");
    }

    [[nodiscard]] std::size_t kernel_size() const noexcept { return kernel_.size(); }

    void reset() {
        std::fill(history_.begin(), history_.end(), T{});
        head_ = 0;
    }

    [[nodiscard]] T process(T sample) {
        history_[head_] = sample;
        T sum{};
        for (std::size_t k = 0; k < kernel_.size(); ++k) {
            const std::size_t index = head_ >= k ? head_ - k : kernel_.size() - (k - head_);
            sum += kernel_[k] * history_[index];
        }
        head_ = head_ + 1 == kernel_.size() ? 0 : head_ + 1;
        return sum;
    }

    [[nodiscard]] std::vector<T> process(std::span<const T> input) {
        std::vector<T> output(input.size());
        for (std::size_t i = 0; i < input.size(); ++i) output[i] = process(input[i]);
        return output;
    }

    [[nodiscard]] std::vector<T> flush() {
        std::vector<T> tail(kernel_.size() - 1);
        for (auto& value : tail) value = process(T{});
        return tail;
    }

private:
    std::vector<T> kernel_;
    std::vector<T> history_;
    std::size_t head_{};
};

template <signal_processing::detail::Sample T>
class UniformPartitionedConvolver {
public:
    UniformPartitionedConvolver(std::span<const T> kernel, std::size_t partition_size)
        : kernel_size_(kernel.size()), partition_size_(partition_size),
          transform_size_(checked_transform_size(kernel_size_, partition_size_)),
          partition_count_(1 + (kernel_size_ - 1) / partition_size_),
          plan_(transform_size_),
          kernel_spectra_(partition_count_, std::vector<Complex>(transform_size_)),
          input_history_(partition_count_, std::vector<Complex>(transform_size_)),
          work_(transform_size_), accumulator_(transform_size_), overlap_(partition_size_, T{}) {
        for (std::size_t p = 0; p < partition_count_; ++p) {
            auto& spectrum = kernel_spectra_[p];
            const std::size_t offset = p * partition_size_;
            const std::size_t count = std::min(partition_size_, kernel_size_ - offset);
            detail::fill_complex<T>(kernel.subspan(offset, count), spectrum);
            plan_.forward_inplace(spectrum);
        }
    }

    [[nodiscard]] std::size_t partition_size() const noexcept { return partition_size_; }
    [[nodiscard]] std::size_t partition_count() const noexcept { return partition_count_; }
    [[nodiscard]] std::size_t transform_size() const noexcept { return transform_size_; }

    void reset() {
        for (auto& spectrum : input_history_)
            std::fill(spectrum.begin(), spectrum.end(), Complex{});
        std::fill(overlap_.begin(), overlap_.end(), T{});
        current_ = 0;
    }

    void process_block(std::span<const T> input, std::span<T> output) {
        if (input.size() != partition_size_ || output.size() != partition_size_)
            throw std::invalid_argument("UniformPartitionedConvolver block size mismatch");

        auto& current_spectrum = input_history_[current_];
        detail::fill_complex<T>(input, current_spectrum);
        plan_.forward_inplace(current_spectrum);
        std::fill(accumulator_.begin(), accumulator_.end(), Complex{});
        for (std::size_t p = 0; p < partition_count_; ++p) {
            const std::size_t history_index =
                current_ >= p ? current_ - p : partition_count_ - (p - current_);
            for (std::size_t k = 0; k < transform_size_; ++k)
                accumulator_[k] += input_history_[history_index][k] * kernel_spectra_[p][k];
        }
        work_ = accumulator_;
        plan_.inverse_inplace(work_);
        for (std::size_t i = 0; i < partition_size_; ++i) {
            output[i] = signal_processing::detail::from_complex<T>(work_[i]) + overlap_[i];
            overlap_[i] = i + partition_size_ < transform_size_
                ? signal_processing::detail::from_complex<T>(work_[i + partition_size_])
                : T{};
        }
        current_ = current_ + 1 == partition_count_ ? 0 : current_ + 1;
    }

private:
    [[nodiscard]] static std::size_t checked_transform_size(std::size_t kernel_size,
                                                             std::size_t partition_size) {
        if (kernel_size == 0)
            throw std::invalid_argument("partitioned convolution requires a nonempty kernel");
        if (partition_size == 0) throw std::invalid_argument("partition size must be nonzero");
        if (partition_size > std::numeric_limits<std::size_t>::max() / 2)
            throw std::length_error("partitioned convolution workspace overflows size_t");
        return detail::next_power_of_two(2 * partition_size);
    }

    using Scalar = signal_processing::detail::scalar_t<T>;
    using Complex = signal_processing::detail::complex_t<T>;
    std::size_t kernel_size_{};
    std::size_t partition_size_{};
    std::size_t transform_size_{};
    std::size_t partition_count_{};
    fft::Radix2Plan<Scalar> plan_;
    std::vector<std::vector<Complex>> kernel_spectra_;
    std::vector<std::vector<Complex>> input_history_;
    std::vector<Complex> work_;
    std::vector<Complex> accumulator_;
    std::vector<T> overlap_;
    std::size_t current_{};
};

// Uniform partitioned frequency-domain convolution. partition_size is the time-domain
// partition length.
template <signal_processing::detail::Sample T>
[[nodiscard]] inline std::vector<T> partitioned(std::span<const T> signal,
                                                 std::span<const T> kernel,
                                                 std::size_t partition_size) {
    const std::size_t output_size = detail::linear_size(signal.size(), kernel.size());
    if (output_size == 0) return {};
    UniformPartitionedConvolver<T> convolver(kernel, partition_size);
    std::vector<T> input_block(partition_size, T{}), output_block(partition_size, T{}), output;
    output.reserve(output_size);
    std::size_t consumed = 0;
    while (output.size() < output_size) {
        std::fill(input_block.begin(), input_block.end(), T{});
        const std::size_t count = std::min(partition_size, signal.size() - std::min(consumed, signal.size()));
        if (count != 0) {
            const auto source = signal.subspan(consumed, count);
            std::copy(source.begin(), source.end(), input_block.begin());
        }
        consumed += count;
        convolver.process_block(input_block, output_block);
        const std::size_t keep = std::min(partition_size, output_size - output.size());
        const auto valid = std::span<const T>(output_block).first(keep);
        output.insert(output.end(), valid.begin(), valid.end());
    }
    return output;
}

}  // namespace signal_processing::convolution
