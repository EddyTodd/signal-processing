// Independent high-precision Fourier reference across float32 and float64,
// including prime/composite sizes, structured/adversarial inputs and inverse
// normalization. No implementation-internal FFT kernels are reused.
#ifdef NDEBUG
#error "numerical_precision must run with C++ assertions enabled in Release"
#endif

#include "signal_processing/fft.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <random>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {
namespace fft = signal_processing::fft;
template <typename T> using C = std::complex<T>;
template <typename T> using Vector = std::vector<C<T>>;

enum class Pattern { random, impulse, dc, alternating, coherent_tone, mixed_dynamic_range };

template <typename T>
Vector<T> stimulus(std::size_t n, Pattern pattern, std::uint32_t seed) {
    Vector<T> values(n);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<T> distribution(T{-1}, T{1});
    for (std::size_t j = 0; j < n; ++j) {
        switch (pattern) {
        case Pattern::random:
            values[j] = {distribution(rng), distribution(rng)};
            break;
        case Pattern::impulse:
            values[j] = j == n/3 ? C<T>{T{1},T{-0.5}} : C<T>{};
            break;
        case Pattern::dc:
            values[j] = {T{0.625},T{-0.125}};
            break;
        case Pattern::alternating:
            values[j] = {j % 2 ? static_cast<T>(-0.8L) : static_cast<T>(0.8L), T{0.0}};
            break;
        case Pattern::coherent_tone: {
            const auto angle = T{2}*std::numbers::pi_v<T>*
                               static_cast<T>((j*3)%n)/static_cast<T>(n);
            values[j] = {std::cos(angle), std::sin(angle)};
            break;
        }
        case Pattern::mixed_dynamic_range:
            values[j] = j % 7 == 0 ? C<T>{T{10},static_cast<T>(0.2L)} :
                C<T>{static_cast<T>(0.0001L)*distribution(rng),static_cast<T>(0.0001L)*distribution(rng)};
            break;
        }
    }
    return values;
}

// This O(N^2) reference accumulates in long double, independently from the
// float/double FFT kernels and from the original binary64 reference suite.
template <typename T>
std::vector<C<long double>> reference_dft(std::span<const C<T>> input,
                                         fft::Direction direction) {
    const std::size_t n = input.size();
    std::vector<C<long double>> result(n);
    const long double sign=direction==fft::Direction::inverse ? 1.L : -1.L;
    for (std::size_t k=0;k<n;++k) {
        C<long double> sum{};
        for (std::size_t j=0;j<n;++j) {
            const long double angle=sign*2.L*std::numbers::pi_v<long double>*
                static_cast<long double>(j)*static_cast<long double>(k)/
                static_cast<long double>(n);
            sum += C<long double>{input[j].real(),input[j].imag()} *
                   C<long double>{std::cos(angle),std::sin(angle)};
        }
        result[k] = direction==fft::Direction::inverse
           ? sum/static_cast<long double>(n) : sum;
    }
    return result;
}

template <typename T>
void compare(std::span<const C<T>> got,
             std::span<const C<long double>> reference,
             const char* algorithm, std::size_t n, Pattern pattern,
             std::uint32_t seed, const char* phase) {
    assert(got.size()==reference.size());
    const long double tol=std::same_as<T,float> ? 0.0008L : 1.0e-10L;
    long double square_error=0.L,square_reference=0.L;
    for(std::size_t i=0;i<got.size();++i){
        const C<long double> actual{got[i].real(),got[i].imag()};
        const auto deviation=std::abs(actual-reference[i]);
        const auto bound=tol*(1.L+std::abs(reference[i]));
        if(!std::isfinite(got[i].real())||!std::isfinite(got[i].imag())||
           !(deviation<=bound)) {
            std::cerr << std::scientific << std::setprecision(17)
                << "PRECISION_FAILURE precision=" << (std::same_as<T,float>?"float32":"float64")
                << " algorithm=" << algorithm << " length=" << n
                << " pattern=" << static_cast<int>(pattern) << " seed=" << seed
                << " stage=" << phase << " index=" << i
                << " got=(" << actual.real() << "," << actual.imag() << ")"
                << " reference=(" << reference[i].real() << "," << reference[i].imag() << ")"
                << " absolute_error=" << deviation << " allowed=" << bound
                << " relative_bound=" << tol << std::endl;
            std::abort();
        }
        square_error += std::norm(actual-reference[i]);
        square_reference += std::norm(reference[i]);
    }
    const auto rms_error=std::sqrt(square_error);
    const auto rms_bound=tol*(1.L+std::sqrt(square_reference));
    if(!(rms_error<=rms_bound)){
        std::cerr << "PRECISION_FAILURE_RMS precision="
            << (std::same_as<T,float>?"float32":"float64")
            << " algorithm=" << algorithm << " length=" << n
            << " pattern=" << static_cast<int>(pattern)
            << " seed=" << seed << " stage=" << phase
            << " rms_error=" << rms_error << " allowed=" << rms_bound << std::endl;
        std::abort();
    }
}

template <typename T>
void check_transform(std::size_t n,Pattern pattern,std::uint32_t seed) {
    const auto input=stimulus<T>(n,pattern,seed);
    const auto reference=reference_dft<T>(input,fft::Direction::forward);
    auto check=[&](const char* name,const Vector<T>& actual){
        compare<T>(actual,reference,name,n,pattern,seed,"forward");
        const auto backward=reference_dft<T>(actual,fft::Direction::inverse);
        compare<T>(input,backward,name,n,pattern,seed,"inverse-reconstruction");
    };
    check("direct",fft::dft<T>(input));
    check("mixed_radix",fft::mixed_radix<T>(input));
    check("bluestein",fft::bluestein<T>(input));
    if(fft::is_power_of_two(n)){
        check("radix2",fft::radix2<T>(input));
        check("radix2_recursive",fft::radix2_recursive<T>(input));
        check("stockham",fft::stockham<T>(input));
        check("radix4",fft::radix4<T>(input));
        check("split_radix",fft::split_radix<T>(input));
        check("modified_split_radix",fft::modified_split_radix<T>(input));
    }
    if(n>=3 && (n==3||n==5||n==7||n==11||n==17||n==31||n==97))
        check("rader",fft::rader<T>(input));
    if(n==6||n==10||n==15||n==21||n==35)
        check("good_thomas",fft::good_thomas<T>(input));

    // Exercise each implementation's inverse normalization independently,
    // instead of checking an inverse that is only computed by the oracle.
    const auto inverse_reference=reference_dft<T>(input,fft::Direction::inverse);
    auto check_inverse=[&](const char* name,const Vector<T>& actual){
        compare<T>(actual,inverse_reference,name,n,pattern,seed,"inverse");
    };
    check_inverse("direct_inverse",fft::dft<T>(input,fft::Direction::inverse));
    check_inverse("mixed_radix_inverse",fft::mixed_radix<T>(input,fft::Direction::inverse));
    check_inverse("bluestein_inverse",fft::bluestein<T>(input,fft::Direction::inverse));
    if(fft::is_power_of_two(n)){
        check_inverse("radix2_inverse",fft::radix2<T>(input,fft::Direction::inverse));
        check_inverse("radix2_recursive_inverse",fft::radix2_recursive<T>(input,fft::Direction::inverse));
        check_inverse("stockham_inverse",fft::stockham<T>(input,fft::Direction::inverse));
        check_inverse("radix4_inverse",fft::radix4<T>(input,fft::Direction::inverse));
        check_inverse("split_radix_inverse",fft::split_radix<T>(input,fft::Direction::inverse));
        check_inverse("modified_split_radix_inverse",fft::modified_split_radix<T>(input,fft::Direction::inverse));
    }
    if(n>=3&&(n==3||n==5||n==7||n==11||n==17||n==31||n==97))
        check_inverse("rader_inverse",fft::rader<T>(input,fft::Direction::inverse));
    if(n==6||n==10||n==15||n==21||n==35)
        check_inverse("good_thomas_inverse",fft::good_thomas<T>(input,0,0,fft::Direction::inverse));
}

template <typename T>
void verify_empty_and_invalid() {
    const Vector<T> empty;
    assert(fft::dft<T>(empty).empty());
    assert(fft::bluestein<T>(empty).empty());
    const Vector<T> nonpower(7,C<T>{1,0});
    bool threw=false;
    try { (void)fft::radix2<T>(nonpower); }
    catch(const std::invalid_argument&) { threw=true; }
    assert(threw);
}

// Independent analytical DC transform identity; the original float32
// precision failure was at length 256, bin 224. Tighten the leakage bound
// separately from the general long-double oracle and include prime lengths.
template <typename T>
void verify_direct_dc_cancellation() {
    constexpr std::array<std::size_t,3> sizes{128U,256U,257U};
    for(const auto n:sizes) {
        const Vector<T> constant(n,C<T>{T{0.625},T{-0.125}});
        for(const auto direction:{fft::Direction::forward,fft::Direction::inverse}) {
            const auto actual=fft::dft<T>(constant,direction);
            const long double scale=direction==fft::Direction::forward
                ? static_cast<long double>(n) : 1.L;
            const C<long double> dc{0.625L*scale,-0.125L*scale};
            const long double tolerance=std::same_as<T,float> ? 5.e-5L : 5.e-8L;
            for(std::size_t k=0;k<n;++k) {
                const C<long double> observed{actual[k].real(),actual[k].imag()};
                const C<long double> expected=k==0 ? dc : C<long double>{};
                const long double error=std::abs(observed-expected);
                const long double bound=tolerance*(1.L+std::abs(expected));
                if(!std::isfinite(actual[k].real())||
                   !std::isfinite(actual[k].imag())||!(error<=bound)) {
                    std::cerr << "PRECISION_FAILURE_DC precision="
                        << (std::same_as<T,float>?"float32":"float64")
                        << " n=" << n << " k=" << k << " direction="
                        << (direction==fft::Direction::forward?"forward":"inverse")
                        << " observed=" << observed << " expected=" << expected
                        << " error=" << error << " bound=" << bound << std::endl;
                    std::abort();
                }
            }
        }
    }
}

template <typename T>
void run_precision_suite() {
    verify_empty_and_invalid<T>();
    verify_direct_dc_cancellation<T>();
    constexpr std::array<std::size_t,24> sizes={1U,2U,3U,4U,5U,6U,7U,8U,10U,
        11U,15U,16U,17U,21U,31U,32U,35U,64U,97U,127U,128U,129U,256U,257U};
    for(const auto n:sizes) {
        for(const auto pattern:{Pattern::random,Pattern::impulse,Pattern::dc,
            Pattern::alternating,Pattern::coherent_tone,Pattern::mixed_dynamic_range})
            check_transform<T>(n,pattern,2718U+static_cast<std::uint32_t>(n));
        // A second independent random stream helps detect seed-specific
        // cancellations without quadratic work on every structured pattern.
        check_transform<T>(n,Pattern::random,998244353U+static_cast<std::uint32_t>(n));
        check_transform<T>(n,Pattern::mixed_dynamic_range,31337U+static_cast<std::uint32_t>(n));
    }
}
} // namespace

int main(){
    run_precision_suite<float>();
    run_precision_suite<double>();
}
