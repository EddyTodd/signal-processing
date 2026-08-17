#pragma once

#include "signal_processing/fft.hpp"

#include <complex>
#include <type_traits>

namespace signal_processing::detail {

template <typename T>
struct SampleTraits;

template <fft::Scalar T>
struct SampleTraits<T> {
    using scalar_type = T;
    static constexpr bool is_complex = false;
};

template <fft::Scalar T>
struct SampleTraits<std::complex<T>> {
    using scalar_type = T;
    static constexpr bool is_complex = true;
};

template <typename T>
concept Sample = requires { typename SampleTraits<T>::scalar_type; };

template <Sample T>
using scalar_t = typename SampleTraits<T>::scalar_type;

template <Sample T>
using complex_t = fft::Complex<scalar_t<T>>;

template <Sample T>
[[nodiscard]] inline complex_t<T> to_complex(const T& value) {
    if constexpr (SampleTraits<T>::is_complex) {
        return value;
    } else {
        return {value, scalar_t<T>{0}};
    }
}

template <Sample T>
[[nodiscard]] inline T from_complex(const complex_t<T>& value) {
    if constexpr (SampleTraits<T>::is_complex) {
        return value;
    } else {
        return value.real();
    }
}

template <Sample T>
[[nodiscard]] inline T conjugate(const T& value) {
    if constexpr (SampleTraits<T>::is_complex) {
        return std::conj(value);
    } else {
        return value;
    }
}

template <Sample T>
[[nodiscard]] inline scalar_t<T> magnitude_squared(const T& value) {
    if constexpr (SampleTraits<T>::is_complex) {
        return std::norm(value);
    } else {
        return value * value;
    }
}

}  // namespace signal_processing::detail
