#pragma once

#include "signal_processing/fft.hpp"
#include "signal_processing/windows.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace signal_processing::fir_design {

// Window-method FIR design.
// https://en.wikipedia.org/wiki/Finite_impulse_response#Window_design_method

namespace detail {

template <fft::Scalar T>
[[nodiscard]] inline T sinc(T x) {
    if (x == T{0}) return T{1};
    const T pix = std::numbers::pi_v<T> * x;
    return std::sin(pix) / pix;
}

template <fft::Scalar T>
void check_frequency(T frequency) {
    if (!(frequency >= T{0} && frequency <= T{0.5}))
        throw std::invalid_argument("normalized frequency must be in [0, 0.5]");
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> apply_window(std::vector<T> coefficients,
                                                  std::span<const T> window) {
    if (window.size() != coefficients.size())
        throw std::invalid_argument("FIR design window length must match tap count");
    for (std::size_t i = 0; i < coefficients.size(); ++i) coefficients[i] *= window[i];
    return coefficients;
}

template <fft::Scalar T>
[[nodiscard]] std::vector<T> solve_linear(std::vector<std::vector<T>> a,
                                           std::vector<T> b) {
    const std::size_t n = b.size();
    if (a.size() != n) throw std::invalid_argument("linear system size mismatch");
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        T best = std::abs(a[col][col]);
        for (std::size_t row = col + 1; row < n; ++row) {
            const T candidate = std::abs(a[row][col]);
            if (candidate > best) { best = candidate; pivot = row; }
        }
        if (best <= T{64} * std::numeric_limits<T>::epsilon())
            throw std::runtime_error("Remez exchange linear system became singular");
        if (pivot != col) { std::swap(a[pivot], a[col]); std::swap(b[pivot], b[col]); }

        const T scale = a[col][col];
        for (std::size_t j = col; j < n; ++j) a[col][j] /= scale;
        b[col] /= scale;
        for (std::size_t row = 0; row < n; ++row) {
            if (row == col) continue;
            const T factor = a[row][col];
            if (factor == T{0}) continue;
            for (std::size_t j = col; j < n; ++j) a[row][j] -= factor * a[col][j];
            b[row] -= factor * b[col];
        }
    }
    return b;
}

template <fft::Scalar T>
struct GridPoint {
    T frequency{};
    T desired{};
    T weight{};
};

template <fft::Scalar T>
[[nodiscard]] std::vector<GridPoint<T>> build_grid(std::span<const T> bands,
                                                    std::span<const T> desired,
                                                    std::span<const T> weights,
                                                    std::size_t density,
                                                    std::size_t cosine_count) {
    if (bands.size() < 2 || (bands.size() & 1U) != 0U)
        throw std::invalid_argument("Remez bands must contain low/high pairs");
    const std::size_t band_count = bands.size() / 2;
    if (desired.size() != band_count || weights.size() != band_count)
        throw std::invalid_argument("Remez desired/weight count must match band count");
    if (density < 4) throw std::invalid_argument("Remez grid density must be at least 4");

    std::vector<GridPoint<T>> grid;
    for (std::size_t band = 0; band < band_count; ++band) {
        const T lo = bands[2 * band];
        const T hi = bands[2 * band + 1];
        check_frequency(lo); check_frequency(hi);
        if (hi < lo || (band != 0 && lo < bands[2 * band - 1]))
            throw std::invalid_argument("Remez bands must be ordered and non-overlapping");
        if (!(weights[band] > T{0}))
            throw std::invalid_argument("Remez weights must be positive");

        const T width = hi - lo;
        const std::size_t points = std::max<std::size_t>(2,
            static_cast<std::size_t>(std::ceil(width * T{2} *
                static_cast<T>(density * cosine_count))) + 1);
        for (std::size_t i = 0; i < points; ++i) {
            if (band != 0 && i == 0 && lo == bands[2 * band - 1]) continue;
            const T fraction = points == 1 ? T{} : static_cast<T>(i) / static_cast<T>(points - 1);
            grid.push_back({lo + width * fraction, desired[band], weights[band]});
        }
    }
    return grid;
}

template <fft::Scalar T>
[[nodiscard]] std::vector<T> exchange_solution(const std::vector<GridPoint<T>>& grid,
                                                std::span<const std::size_t> extrema,
                                                std::size_t cosine_count) {
    const std::size_t equations = cosine_count + 1;
    std::vector<std::vector<T>> matrix(equations, std::vector<T>(equations));
    std::vector<T> rhs(equations);
    for (std::size_t i = 0; i < equations; ++i) {
        const auto& point = grid[extrema[i]];
        for (std::size_t k = 0; k < cosine_count; ++k)
            matrix[i][k] = std::cos(T{2} * std::numbers::pi_v<T> * point.frequency *
                                    static_cast<T>(k));
        matrix[i][cosine_count] = ((i & 1U) == 0U ? T{1} : T{-1}) / point.weight;
        rhs[i] = point.desired;
    }
    return solve_linear<T>(std::move(matrix), std::move(rhs));
}

template <fft::Scalar T>
[[nodiscard]] T weighted_error(const GridPoint<T>& point, std::span<const T> solution,
                                std::size_t cosine_count) {
    T response{};
    for (std::size_t k = 0; k < cosine_count; ++k)
        response += solution[k] * std::cos(T{2} * std::numbers::pi_v<T> * point.frequency *
                                           static_cast<T>(k));
    return point.weight * (response - point.desired);
}

template <fft::Scalar T>
[[nodiscard]] std::vector<std::size_t> select_extrema(const std::vector<T>& error,
                                                       std::size_t count) {
    std::vector<std::size_t> candidates;
    if (error.empty()) return candidates;
    candidates.push_back(0);
    for (std::size_t i = 1; i + 1 < error.size(); ++i) {
        const T magnitude = std::abs(error[i]);
        if (magnitude >= std::abs(error[i - 1]) && magnitude >= std::abs(error[i + 1]))
            candidates.push_back(i);
    }
    if (error.size() > 1) candidates.push_back(error.size() - 1);

    std::vector<std::size_t> alternating;
    for (const auto index : candidates) {
        if (alternating.empty()) { alternating.push_back(index); continue; }
        const bool same_sign = std::signbit(error[alternating.back()]) == std::signbit(error[index]);
        if (same_sign) {
            if (std::abs(error[index]) > std::abs(error[alternating.back()])) alternating.back() = index;
        } else {
            alternating.push_back(index);
        }
    }

    while (alternating.size() > count) {
        auto it = std::min_element(alternating.begin(), alternating.end(), [&](auto a, auto b) {
            return std::abs(error[a]) < std::abs(error[b]);
        });
        alternating.erase(it);
        std::vector<std::size_t> collapsed;
        for (const auto index : alternating) {
            if (!collapsed.empty() &&
                std::signbit(error[collapsed.back()]) == std::signbit(error[index])) {
                if (std::abs(error[index]) > std::abs(error[collapsed.back()])) collapsed.back() = index;
            } else collapsed.push_back(index);
        }
        alternating = std::move(collapsed);
    }

    if (alternating.size() < count) {
        std::vector<std::size_t> order(error.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](auto a, auto b) {
            return std::abs(error[a]) > std::abs(error[b]);
        });
        for (const auto index : order) {
            if (std::find(alternating.begin(), alternating.end(), index) != alternating.end()) continue;
            alternating.push_back(index);
            if (alternating.size() == count) break;
        }
        std::sort(alternating.begin(), alternating.end());
    }
    return alternating;
}

}  // namespace detail

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> lowpass_windowed_sinc(std::size_t tap_count, T cutoff,
                                                           std::span<const T> window) {
    if (tap_count == 0) throw std::invalid_argument("FIR design requires at least one tap");
    detail::check_frequency(cutoff);
    if (!(cutoff > T{0} && cutoff < T{0.5}))
        throw std::invalid_argument("low-pass cutoff must lie strictly inside (0, 0.5)");
    const T center = static_cast<T>(tap_count - 1) / T{2};
    std::vector<T> coefficients(tap_count);
    for (std::size_t n = 0; n < tap_count; ++n) {
        const T offset = static_cast<T>(n) - center;
        coefficients[n] = T{2} * cutoff * detail::sinc<T>(T{2} * cutoff * offset);
    }
    return detail::apply_window<T>(std::move(coefficients), window);
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> lowpass_windowed_sinc(std::size_t tap_count, T cutoff) {
    const auto window = windows::hann<T>(tap_count);
    return lowpass_windowed_sinc<T>(tap_count, cutoff, window);
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> highpass_windowed_sinc(std::size_t tap_count, T cutoff,
                                                            std::span<const T> window) {
    if ((tap_count & 1U) == 0U)
        throw std::invalid_argument("spectral-inversion high-pass design requires an odd tap count");
    auto coefficients = lowpass_windowed_sinc<T>(tap_count, cutoff, window);
    for (auto& value : coefficients) value = -value;
    coefficients[tap_count / 2] += T{1};
    return coefficients;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> bandpass_windowed_sinc(std::size_t tap_count, T low_cutoff,
                                                            T high_cutoff,
                                                            std::span<const T> window) {
    detail::check_frequency(low_cutoff); detail::check_frequency(high_cutoff);
    if (!(T{0} < low_cutoff && low_cutoff < high_cutoff && high_cutoff < T{0.5}))
        throw std::invalid_argument("band-pass cutoffs must satisfy 0 < low < high < 0.5");
    auto high = lowpass_windowed_sinc<T>(tap_count, high_cutoff, window);
    const auto low = lowpass_windowed_sinc<T>(tap_count, low_cutoff, window);
    for (std::size_t i = 0; i < tap_count; ++i) high[i] -= low[i];
    return high;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> bandstop_windowed_sinc(std::size_t tap_count, T low_cutoff,
                                                            T high_cutoff,
                                                            std::span<const T> window) {
    if ((tap_count & 1U) == 0U)
        throw std::invalid_argument("spectral-inversion band-stop design requires an odd tap count");
    auto coefficients = bandpass_windowed_sinc<T>(tap_count, low_cutoff, high_cutoff, window);
    for (auto& value : coefficients) value = -value;
    coefficients[tap_count / 2] += T{1};
    return coefficients;
}

// Type-I odd-length Parks-McClellan / Remez exchange design.
// https://en.wikipedia.org/wiki/Remez_algorithm
// https://en.wikipedia.org/wiki/Parks%E2%80%93McClellan_filter_design_algorithm

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> remez_type1(std::size_t tap_count,
                                                 std::span<const T> bands,
                                                 std::span<const T> desired,
                                                 std::span<const T> weights,
                                                 std::size_t grid_density = 16,
                                                 std::size_t max_iterations = 64) {
    if (tap_count < 3 || (tap_count & 1U) == 0U)
        throw std::invalid_argument("remez_type1 requires an odd tap count >= 3");
    const std::size_t cosine_count = (tap_count + 1) / 2;
    const std::size_t extremum_count = cosine_count + 1;
    const auto grid = detail::build_grid<T>(bands, desired, weights, grid_density, cosine_count);
    if (grid.size() < extremum_count)
        throw std::invalid_argument("Remez grid is too small for the requested tap count");

    std::vector<std::size_t> extrema(extremum_count);
    for (std::size_t i = 0; i < extremum_count; ++i)
        extrema[i] = i * (grid.size() - 1) / (extremum_count - 1);

    std::vector<T> solution;
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        solution = detail::exchange_solution<T>(grid, extrema, cosine_count);
        std::vector<T> error(grid.size());
        for (std::size_t i = 0; i < grid.size(); ++i)
            error[i] = detail::weighted_error<T>(grid[i], solution, cosine_count);
        auto next = detail::select_extrema<T>(error, extremum_count);
        if (next.size() != extremum_count)
            throw std::runtime_error("Remez exchange failed to find enough extrema");
        if (next == extrema) break;
        extrema = std::move(next);
    }
    solution = detail::exchange_solution<T>(grid, extrema, cosine_count);

    std::vector<T> coefficients(tap_count, T{});
    const std::size_t center = tap_count / 2;
    coefficients[center] = solution[0];
    for (std::size_t k = 1; k < cosine_count; ++k) {
        const T value = solution[k] / T{2};
        coefficients[center - k] = value;
        coefficients[center + k] = value;
    }
    return coefficients;
}

template <fft::Scalar T>
[[nodiscard]] inline std::vector<T> remez_lowpass(std::size_t tap_count, T pass_edge,
                                                   T stop_edge, T pass_weight = T{1},
                                                   T stop_weight = T{1},
                                                   std::size_t grid_density = 16) {
    if (!(T{0} <= pass_edge && pass_edge < stop_edge && stop_edge <= T{0.5}))
        throw std::invalid_argument("low-pass edges must satisfy 0 <= pass < stop <= 0.5");
    const std::vector<T> bands{T{0}, pass_edge, stop_edge, T{0.5}};
    const std::vector<T> desired{T{1}, T{0}};
    const std::vector<T> weights{pass_weight, stop_weight};
    return remez_type1<T>(tap_count, bands, desired, weights, grid_density);
}

}  // namespace signal_processing::fir_design
