// Independent textbook long-double Goertzel and Hilbert/analytic-signal oracles.
// Unlike the production algorithms, these references use direct O(N^2) DFTs.
#ifdef NDEBUG
#error "Analytic signal oracles require active Release assertions"
#endif
#include "signal_processing/goertzel.hpp"
#include "signal_processing/hilbert.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
using Wide=long double;
using CW=std::complex<Wide>;
constexpr Wide tau=Wide{2}*std::numbers::pi_v<Wide>;

template <typename Sample>
using Scalar=signal_processing::detail::scalar_t<Sample>;

template <typename Sample>
CW widen(Sample v) {
    if constexpr (std::same_as<Sample,Scalar<Sample>>)
        return {static_cast<Wide>(v), Wide{}};
    else
        return {static_cast<Wide>(v.real()),static_cast<Wide>(v.imag())};
}

template <typename Sample>
std::vector<Sample> make_input(std::size_t n,unsigned pattern) {
    std::vector<Sample> out(n);
    using T=Scalar<Sample>;
    for(std::size_t i=0;i<n;++i) {
        const Wide j=static_cast<Wide>(i);
        const Wide re=pattern==0 ? Wide{0.37}*std::cos(Wide{0.21}*j)-
                                    Wide{0.19}*std::sin(Wide{0.74}*j) :
                        pattern==1 ? (i==0?Wide{1}:Wide{}) :
                                     static_cast<Wide>(static_cast<int>((7*i+3)%13)-6)/Wide{9};
        if constexpr (std::same_as<Sample,T>)
            out[i]=static_cast<T>(re);
        else
            out[i]=Sample{static_cast<T>(re),
                          static_cast<T>(pattern==1?Wide{}:
                            Wide{0.16}*std::sin(Wide{0.29}*j)+Wide{0.07})};
    }
    return out;
}

void check(CW actual,CW expected,Wide eps,
           const char* method,std::size_t n,Wide frequency,
           std::size_t index) {
    const Wide error=std::abs(actual-expected);
    const Wide bound=eps*(Wide{1}+std::abs(expected));
    if(!std::isfinite(actual.real())||!std::isfinite(actual.imag())||
       !(error<=bound)) {
        std::cerr<<"ANALYTIC_ORACLE_FAILURE method="<<method<<" n="<<n
                 <<" f="<<frequency<<" index="<<index
                 <<" actual="<<actual<<" expected="<<expected
                 <<" error="<<error<<" bound="<<bound<<'\n';
        std::abort();
    }
}

template <typename Sample>
CW direct_frequency(const std::vector<Sample>& input,Wide frequency) {
    CW result{};
    for(std::size_t i=0;i<input.size();++i) {
        const Wide angle=-tau*frequency*static_cast<Wide>(i);
        result+=widen(input[i])*CW{std::cos(angle),std::sin(angle)};
    }
    return result;
}

template <typename Sample>
void check_goertzel() {
    using T=Scalar<Sample>;
    const Wide tolerance=std::same_as<T,float>?Wide{0.008}:Wide{3.e-9};
    for(const std::size_t n:{0U,1U,2U,3U,5U,8U,11U,16U,31U,47U,64U,97U})
        for(unsigned pattern=0;pattern<3;++pattern) {
            const auto input=make_input<Sample>(n,pattern);
            for(const Wide freq:{Wide{},Wide{0.017},Wide{0.137},
                                 Wide{0.31},Wide{0.499},Wide{0.75},
                                 Wide{0.98},Wide{0.999}}) {
                const T f=static_cast<T>(freq);
                check(widen(signal_processing::goertzel::frequency<Sample>(input,f)),
                      direct_frequency(input,static_cast<Wide>(f)),
                      tolerance,"off_bin",n,static_cast<Wide>(f),pattern);
            }
            if(n==0) {
                assert(signal_processing::goertzel::bin<Sample>(input,0)==
                       std::complex<T>{});
            } else {
                for(std::size_t k=0;k<n;++k) {
                    const Wide f=static_cast<Wide>(k)/static_cast<Wide>(n);
                    check(widen(signal_processing::goertzel::bin<Sample>(input,k)),
                          direct_frequency(input,f),
                          tolerance,"dft_bin",n,f,k);
                }
            }
        }
    const std::vector<Sample> input=make_input<Sample>(5,0);
    for(const T bad:{std::numeric_limits<T>::quiet_NaN(),
                     std::numeric_limits<T>::infinity(),
                     -std::numeric_limits<T>::infinity(),T{-1},T{1}}) {
        bool threw=false;
        try {(void)signal_processing::goertzel::frequency<Sample>(input,bad);}
        catch(const std::invalid_argument&) {threw=true;}
        assert(threw);
    }
    bool threw=false;
    try {(void)signal_processing::goertzel::bin<Sample>(input,input.size());}
    catch(const std::out_of_range&) {threw=true;}
    assert(threw);
    // Nonfinite real or imaginary input values must not escape from
    // the real/complex Goertzel API as NaN/Inf coefficients.
    for (const T bad:{std::numeric_limits<T>::quiet_NaN(),
                     std::numeric_limits<T>::infinity(),
                     -std::numeric_limits<T>::infinity()}) {
        for (unsigned component=0;component<
             (std::same_as<Sample,T>?1U:2U);++component) {
            auto corrupted=input;
            if constexpr (std::same_as<Sample,T>) {
                corrupted[1]=bad;
            } else {
                corrupted[1]=component==0?Sample{bad,T{}}:Sample{T{},bad};
            }
            bool input_rejected=false;
            try {(void)signal_processing::goertzel::frequency<Sample>(
                     corrupted,T{0.23});}
            catch(const std::invalid_argument&) {input_rejected=true;}
            assert(input_rejected);
            input_rejected=false;
            try {(void)signal_processing::goertzel::bin<Sample>(
                     corrupted,1);}
            catch(const std::invalid_argument&) {input_rejected=true;}
            assert(input_rejected);
        }
    }
    // A representable input can still blow up the unnormalized
    // recurrence at zero frequency. Reject instead of returning Inf.
    const T huge=std::numeric_limits<T>::max()/T{2};
    const std::vector<Sample> huge_input(4,Sample{huge});
    bool overflow_rejected=false;
    try {(void)signal_processing::goertzel::frequency<Sample>(
             huge_input,T{});}
    catch(const std::overflow_error&) {overflow_rejected=true;}
    assert(overflow_rejected);
}

template <typename T>
std::vector<CW> direct_analytic(const std::vector<T>& input) {
    const std::size_t n=input.size();
    if(n==0)return {};
    std::vector<CW> spectrum(n,CW{});
    for(std::size_t k=0;k<n;++k) {
        for(std::size_t t=0;t<n;++t) {
            const Wide a=-tau*static_cast<Wide>(k)*static_cast<Wide>(t)/
                          static_cast<Wide>(n);
            spectrum[k]+=static_cast<Wide>(input[t])*
                         CW{std::cos(a),std::sin(a)};
        }
        const bool positive=k>0&&k<(n+1)/2;
        const bool negative=k>n/2;
        if(positive)spectrum[k]*=Wide{2};
        if(negative)spectrum[k]=CW{};
    }
    std::vector<CW> expected(n,CW{});
    for(std::size_t t=0;t<n;++t) {
        for(std::size_t k=0;k<n;++k) {
            const Wide a=tau*static_cast<Wide>(k)*static_cast<Wide>(t)/
                         static_cast<Wide>(n);
            expected[t]+=spectrum[k]*CW{std::cos(a),std::sin(a)};
        }
        expected[t]/=static_cast<Wide>(n);
    }
    return expected;
}

template <typename T>
void check_hilbert() {
    const Wide tolerance=std::same_as<T,float>?Wide{0.0009}:Wide{3.e-10};
    for(const std::size_t n:{0U,1U,2U,3U,4U,5U,7U,8U,9U,16U,17U,32U,63U})
        for(unsigned pattern=0;pattern<3;++pattern) {
            const auto input=make_input<T>(n,pattern);
            const auto expected=direct_analytic(input);
            const auto analytic=signal_processing::hilbert::analytic_bluestein<T>(input);
            const auto imag=signal_processing::hilbert::transform_bluestein<T>(input);
            const auto env=signal_processing::hilbert::envelope_bluestein<T>(input);
            assert(analytic.size()==n&&imag.size()==n&&env.size()==n);
            for(std::size_t i=0;i<n;++i) {
                check(widen(analytic[i]),expected[i],tolerance,
                      "analytic",n,Wide{},i);
                check(CW{static_cast<Wide>(imag[i]),Wide{}},
                      CW{expected[i].imag(),Wide{}},tolerance,
                      "hilbert_imag",n,Wide{},i);
                check(CW{static_cast<Wide>(env[i]),Wide{}},
                      CW{std::abs(expected[i]),Wide{}},tolerance,
                      "envelope",n,Wide{},i);
                check(CW{static_cast<Wide>(analytic[i].real()),Wide{}},
                      CW{static_cast<Wide>(input[i]),Wide{}},tolerance,
                      "real_reconstruction",n,Wide{},i);
            }
        }
    for (const T bad:{std::numeric_limits<T>::quiet_NaN(),
                     std::numeric_limits<T>::infinity(),
                     -std::numeric_limits<T>::infinity()}) {
        auto invalid=make_input<T>(7,0);
        invalid[3]=bad;
        bool rejected=false;
        try {(void)signal_processing::hilbert::analytic_bluestein<T>(invalid);}
        catch(const std::invalid_argument&) {rejected=true;}
        assert(rejected);
        rejected=false;
        try {(void)signal_processing::hilbert::transform_bluestein<T>(invalid);}
        catch(const std::invalid_argument&) {rejected=true;}
        assert(rejected);
        rejected=false;
        try {(void)signal_processing::hilbert::envelope_bluestein<T>(invalid);}
        catch(const std::invalid_argument&) {rejected=true;}
        assert(rejected);
    }
}
}
int main() {
    check_goertzel<float>();
    check_goertzel<double>();
    check_goertzel<std::complex<float>>();
    check_goertzel<std::complex<double>>();
    check_hilbert<float>();
    check_hilbert<double>();
}
