#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#if (defined(__GNUC__) || defined(__clang__)) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#define SIGNAL_PROCESSING_FFT_X86_TARGETS 1
#else
#define SIGNAL_PROCESSING_FFT_X86_TARGETS 0
#endif

namespace signal_processing::fft {

enum class KernelIsa { scalar, avx2, avx512 };

struct KernelCapabilities {
    bool avx2{};
    bool avx512{};
};

namespace detail::simd {
using Complex64 = Complex<double>;
using Swap = std::pair<std::size_t, std::size_t>;
using KernelFn = void (*)(Complex64*, std::size_t, const Swap*, std::size_t,
                          const std::size_t*, std::size_t, const Complex64*, bool);

inline void apply_permutation(Complex64* data, const Swap* swaps, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) std::swap(data[swaps[i].first], data[swaps[i].second]);
}

inline void scalar_kernel(Complex64* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                          const std::size_t* offsets, std::size_t stage_count,
                          const Complex64* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);
    std::size_t stage = 0;
    for (std::size_t length = 2; stage < stage_count; ++stage, length <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = length / 2;
        for (std::size_t base = 0; base < n; base += length) {
            for (std::size_t offset = 0; offset < half; ++offset) {
                const Complex64 w = inverse ? std::conj(stage_twiddles[offset]) : stage_twiddles[offset];
                const auto upper = data[base + offset];
                const auto lower = data[base + offset + half] * w;
                data[base + offset] = upper + lower;
                data[base + offset + half] = upper - lower;
            }
        }
    }
    if (inverse) {
        const double scale = 1.0 / static_cast<double>(n);
        for (std::size_t i = 0; i < n; ++i) data[i] *= scale;
    }
}

#if SIGNAL_PROCESSING_FFT_X86_TARGETS
__attribute__((target("avx2,fma")))
inline void avx2_kernel(Complex64* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                        const std::size_t* offsets, std::size_t stage_count,
                        const Complex64* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);
    const __m256d conjugate_mask = _mm256_setr_pd(1.0, -1.0, 1.0, -1.0);
    std::size_t stage = 0;
    for (std::size_t length = 2; stage < stage_count; ++stage, length <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = length / 2;
        for (std::size_t base = 0; base < n; base += length) {
            std::size_t offset = 0;
            for (; offset + 2 <= half; offset += 2) {
                const auto upper = _mm256_loadu_pd(reinterpret_cast<const double*>(data + base + offset));
                const auto lower = _mm256_loadu_pd(reinterpret_cast<const double*>(data + base + offset + half));
                auto w = _mm256_loadu_pd(reinterpret_cast<const double*>(stage_twiddles + offset));
                if (inverse) w = _mm256_mul_pd(w, conjugate_mask);
                const auto wr = _mm256_movedup_pd(w);
                const auto wi = _mm256_permute_pd(w, 0xF);
                const auto swapped = _mm256_permute_pd(lower, 0x5);
                const auto cross = _mm256_mul_pd(swapped, wi);
                const auto product = _mm256_fmaddsub_pd(lower, wr, cross);
                _mm256_storeu_pd(reinterpret_cast<double*>(data + base + offset), _mm256_add_pd(upper, product));
                _mm256_storeu_pd(reinterpret_cast<double*>(data + base + offset + half), _mm256_sub_pd(upper, product));
            }
            for (; offset < half; ++offset) {
                const Complex64 w = inverse ? std::conj(stage_twiddles[offset]) : stage_twiddles[offset];
                const auto upper = data[base + offset];
                const auto lower = data[base + offset + half] * w;
                data[base + offset] = upper + lower;
                data[base + offset + half] = upper - lower;
            }
        }
    }
    if (inverse) {
        const auto scale = _mm256_set1_pd(1.0 / static_cast<double>(n));
        std::size_t i = 0;
        for (; i + 2 <= n; i += 2) {
            const auto value = _mm256_loadu_pd(reinterpret_cast<const double*>(data + i));
            _mm256_storeu_pd(reinterpret_cast<double*>(data + i), _mm256_mul_pd(value, scale));
        }
        for (; i < n; ++i) data[i] /= static_cast<double>(n);
    }
}

__attribute__((target("avx512f,avx512dq,avx512vl,fma")))
inline void avx512_kernel(Complex64* data, std::size_t n, const Swap* swaps, std::size_t swap_count,
                          const std::size_t* offsets, std::size_t stage_count,
                          const Complex64* twiddles, bool inverse) {
    apply_permutation(data, swaps, swap_count);
    const __m512d multiply_sign = _mm512_setr_pd(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0);
    const __m512d conjugate_mask = _mm512_setr_pd(1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0);
    std::size_t stage = 0;
    for (std::size_t length = 2; stage < stage_count; ++stage, length <<= 1) {
        const auto* stage_twiddles = twiddles + offsets[stage];
        const auto half = length / 2;
        for (std::size_t base = 0; base < n; base += length) {
            std::size_t offset = 0;
            for (; offset + 4 <= half; offset += 4) {
                const auto upper = _mm512_loadu_pd(reinterpret_cast<const void*>(data + base + offset));
                const auto lower = _mm512_loadu_pd(reinterpret_cast<const void*>(data + base + offset + half));
                auto w = _mm512_loadu_pd(reinterpret_cast<const void*>(stage_twiddles + offset));
                if (inverse) w = _mm512_mul_pd(w, conjugate_mask);
                const auto wr = _mm512_movedup_pd(w);
                const auto wi = _mm512_permute_pd(w, 0xFF);
                const auto swapped = _mm512_permute_pd(lower, 0x55);
                const auto cross = _mm512_mul_pd(_mm512_mul_pd(swapped, wi), multiply_sign);
                const auto product = _mm512_fmadd_pd(lower, wr, cross);
                _mm512_storeu_pd(reinterpret_cast<void*>(data + base + offset), _mm512_add_pd(upper, product));
                _mm512_storeu_pd(reinterpret_cast<void*>(data + base + offset + half), _mm512_sub_pd(upper, product));
            }
            for (; offset < half; ++offset) {
                const Complex64 w = inverse ? std::conj(stage_twiddles[offset]) : stage_twiddles[offset];
                const auto upper = data[base + offset];
                const auto lower = data[base + offset + half] * w;
                data[base + offset] = upper + lower;
                data[base + offset + half] = upper - lower;
            }
        }
    }
    if (inverse) {
        const auto scale = _mm512_set1_pd(1.0 / static_cast<double>(n));
        std::size_t i = 0;
        for (; i + 4 <= n; i += 4) {
            const auto value = _mm512_loadu_pd(reinterpret_cast<const void*>(data + i));
            _mm512_storeu_pd(reinterpret_cast<void*>(data + i), _mm512_mul_pd(value, scale));
        }
        for (; i < n; ++i) data[i] /= static_cast<double>(n);
    }
}
#endif

[[nodiscard]] inline KernelCapabilities detect_capabilities() noexcept {
    KernelCapabilities caps{};
#if SIGNAL_PROCESSING_FFT_X86_TARGETS
    __builtin_cpu_init();
    caps.avx2 = __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
    caps.avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
                  __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("fma");
#endif
    return caps;
}

[[nodiscard]] inline KernelFn function_for(KernelIsa isa) {
    if (isa == KernelIsa::scalar) return scalar_kernel;
#if SIGNAL_PROCESSING_FFT_X86_TARGETS
    if (isa == KernelIsa::avx2) return avx2_kernel;
    if (isa == KernelIsa::avx512) return avx512_kernel;
#endif
    throw std::invalid_argument("requested FFT SIMD kernel is not compiled for this target");
}
}  // namespace detail::simd

[[nodiscard]] inline KernelCapabilities kernel_capabilities() noexcept {
    static const KernelCapabilities capabilities = detail::simd::detect_capabilities();
    return capabilities;
}

template <KernelIsa Isa>
[[nodiscard]] inline bool kernel_available() noexcept {
    if constexpr (Isa == KernelIsa::scalar) return true;
    const auto caps = kernel_capabilities();
    if constexpr (Isa == KernelIsa::avx2) return caps.avx2;
    return caps.avx512;
}

class KernelRadix2Plan {
public:
    explicit KernelRadix2Plan(std::size_t n, KernelIsa isa = KernelIsa::scalar)
        : n_(n), selected_(isa) {
        if (n == 0 || !detail::is_power_of_two(n))
            throw std::invalid_argument("KernelRadix2Plan requires a power-of-two length >= 1");
        const auto bits = detail::ilog2(n);
        swaps_.reserve(n / 2);
        for (std::size_t i = 0; i < n; ++i) {
            auto value = i; std::size_t reversed = 0;
            for (std::size_t bit = 0; bit < bits; ++bit) {
                reversed = (reversed << 1) | (value & 1U); value >>= 1;
            }
            if (i < reversed) swaps_.emplace_back(i, reversed);
        }
        stage_offsets_.reserve(bits);
        std::size_t total = 0;
        for (std::size_t length = 2; length <= n;) {
            stage_offsets_.push_back(total); total += length / 2;
            if (length == n) break;
            length <<= 1;
        }
        twiddles_.resize(total);
        for (std::size_t stage = 0, length = 2; stage < stage_offsets_.size(); ++stage, length <<= 1) {
            for (std::size_t offset = 0; offset < length / 2; ++offset) {
                const double angle = -2.0 * std::numbers::pi_v<double> *
                    static_cast<double>(offset) / static_cast<double>(length);
                twiddles_[stage_offsets_[stage] + offset] = {std::cos(angle), std::sin(angle)};
            }
        }
        const auto caps = kernel_capabilities();
        if (isa == KernelIsa::avx2 && !caps.avx2) throw std::invalid_argument("AVX2/FMA FFT kernel unavailable");
        if (isa == KernelIsa::avx512 && !caps.avx512) throw std::invalid_argument("AVX-512/FMA FFT kernel unavailable");
        fn_ = detail::simd::function_for(isa);
    }

    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] KernelIsa selected_isa() const noexcept { return selected_; }
    [[nodiscard]] std::size_t scratch_size() const noexcept { return 0; }
    [[nodiscard]] std::size_t stored_twiddles() const noexcept { return twiddles_.size(); }
    [[nodiscard]] std::size_t stored_swaps() const noexcept { return swaps_.size(); }

    void forward_inplace(std::span<Complex<double>> data) const { execute(data, false); }
    void inverse_inplace(std::span<Complex<double>> data) const { execute(data, true); }

private:
    void execute(std::span<Complex<double>> data, bool inverse) const {
        if (data.size() != n_) throw std::invalid_argument("KernelRadix2Plan buffer size mismatch");
        if (n_ <= 1) return;
        fn_(data.data(), n_, swaps_.data(), swaps_.size(), stage_offsets_.data(),
            stage_offsets_.size(), twiddles_.data(), inverse);
    }
    detail::simd::KernelFn fn_{};
    std::size_t n_{};
    KernelIsa selected_{KernelIsa::scalar};
    std::vector<detail::simd::Swap> swaps_;
    std::vector<std::size_t> stage_offsets_;
    Vector<double> twiddles_;
};

}  // namespace signal_processing::fft

#undef SIGNAL_PROCESSING_FFT_X86_TARGETS
