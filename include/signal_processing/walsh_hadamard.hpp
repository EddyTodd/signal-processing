#pragma once

#include <bit>
#include <concepts>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace signal_processing::walsh_hadamard {

// Sylvester-ordered Walsh-Hadamard transform.
// https://en.wikipedia.org/wiki/Fast_Walsh%E2%80%93Hadamard_transform

template <typename T>
concept Scalar = std::same_as<T, float> || std::same_as<T, double>;

[[nodiscard]] constexpr bool supported_size(std::size_t n) noexcept {
    return n != 0 && (n & (n - 1)) == 0;
}

template <Scalar T>
[[nodiscard]] inline std::vector<T> direct(std::span<const T> input) {
    const std::size_t n = input.size();
    if (n == 0) return {};
    if (!supported_size(n)) {
        throw std::invalid_argument("Walsh-Hadamard transform requires a power-of-two length");
    }

    using UnsignedSize = std::make_unsigned_t<std::size_t>;
    std::vector<T> output(n, T{});
    for (std::size_t k = 0; k < n; ++k) {
        T sum{};
        for (std::size_t j = 0; j < n; ++j) {
            const auto parity = std::popcount(static_cast<UnsignedSize>(j & k)) & 1;
            sum += parity == 0 ? input[j] : -input[j];
        }
        output[k] = sum;
    }
    return output;
}

template <Scalar T>
inline void fast_inplace(std::span<T> data) {
    const std::size_t n = data.size();
    if (n == 0) return;
    if (!supported_size(n)) {
        throw std::invalid_argument("fast Walsh-Hadamard transform requires a power-of-two length");
    }

    for (std::size_t width = 1; width < n; width <<= 1) {
        const std::size_t block = width << 1;
        for (std::size_t base = 0; base < n; base += block) {
            for (std::size_t j = 0; j < width; ++j) {
                const T a = data[base + j];
                const T b = data[base + j + width];
                data[base + j] = a + b;
                data[base + j + width] = a - b;
            }
        }
    }
}

template <Scalar T>
[[nodiscard]] inline std::vector<T> fast(std::span<const T> input) {
    std::vector<T> output(input.begin(), input.end());
    fast_inplace<T>(output);
    return output;
}

}  // namespace signal_processing::walsh_hadamard
