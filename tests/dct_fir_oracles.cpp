// Independent DCT-I/II/III/IV and FIR textbook oracles.
#ifdef NDEBUG
#error "DCT/FIR oracles must run with Release assertions enabled"
#endif

#include "signal_processing/dct.hpp"
#include "signal_processing/fir.hpp"

#include <algorithm>
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
#include <span>
#include <stdexcept>
#include <vector>

namespace {
namespace dct = signal_processing::dct;
namespace fir = signal_processing::fir;
using Wide = long double;
using ComplexWide = std::complex<Wide>;

template <typename Sample>
ComplexWide widen(Sample x) {
    if constexpr (std::same_as<Sample,float> || std::same_as<Sample,double>)
        return {static_cast<Wide>(x),Wide{}};
    else
        return {static_cast<Wide>(x.real()),static_cast<Wide>(x.imag())};
}

void check(ComplexWide actual,ComplexWide expected,Wide tolerance,
           const char* method,std::size_t n,std::size_t m,std::size_t pattern,
           std::size_t index) {
    const Wide error=std::abs(actual-expected);
    const Wide bound=tolerance*(Wide{1}+std::abs(expected));
    if(!std::isfinite(actual.real())||!std::isfinite(actual.imag())||
       !(error<=bound)) {
        std::cerr << "DCT_FIR_ORACLE_FAILURE method=" << method
                  << " n=" << n << " m=" << m
                  << " pattern=" << pattern << " index=" << index
                  << " actual=" << actual << " expected=" << expected
                  << " error=" << error << " bound=" << bound << '\n';
        std::abort();
    }
}

template <typename T>
std::vector<T> signal(std::size_t n,std::size_t pattern) {
    std::vector<T> x(n);
    for(std::size_t j=0;j<n;++j) {
        Wide value{};
        if(pattern==0) value=Wide{1};
        else if(pattern==1) value=j==0?Wide{1}:Wide{};
        else if(pattern==2) value=j%2?Wide{-1}:Wide{1};
        else if(pattern==3) value=std::sin(Wide{0.19}*static_cast<Wide>(j))+
                                   Wide{0.13}*std::cos(Wide{0.71}*static_cast<Wide>(j));
        else value=Wide{0.03125}*static_cast<Wide>(static_cast<int>(j%7)-3);
        x[j]=static_cast<T>(value);
    }
    return x;
}

// Defining unnormalized textbook cosine transforms; all computations use
// long double rather than library DCT or FFT kernels.
template <typename T>
std::vector<Wide> reference_dct(std::span<const T> x,int kind) {
    const std::size_t n=x.size();
    std::vector<Wide> out(n,Wide{});
    if(n==0)return out;
    if(kind==1 && n==1) {
        out[0]=static_cast<Wide>(x[0]);
        return out;
    }
    const Wide pi=std::numbers::pi_v<Wide>;
    for(std::size_t k=0;k<n;++k) {
        Wide sum{};
        if(kind==1) {
            sum=static_cast<Wide>(x[0])+
                (k%2==0?Wide{1}:Wide{-1})*static_cast<Wide>(x[n-1]);
            for(std::size_t j=1;j+1<n;++j)
                sum+=Wide{2}*static_cast<Wide>(x[j])*
                     std::cos(pi*static_cast<Wide>(j)*static_cast<Wide>(k)/
                              static_cast<Wide>(n-1));
        } else if(kind==2) {
            for(std::size_t j=0;j<n;++j)
                sum+=Wide{2}*static_cast<Wide>(x[j])*
                     std::cos(pi*(static_cast<Wide>(j)+Wide{0.5})*
                              static_cast<Wide>(k)/static_cast<Wide>(n));
        } else if(kind==3) {
            sum=static_cast<Wide>(x[0]);
            for(std::size_t j=1;j<n;++j)
                sum+=Wide{2}*static_cast<Wide>(x[j])*
                     std::cos(pi*static_cast<Wide>(j)*
                              (static_cast<Wide>(k)+Wide{0.5})/
                              static_cast<Wide>(n));
        } else {
            for(std::size_t j=0;j<n;++j)
                sum+=Wide{2}*static_cast<Wide>(x[j])*
                     std::cos(pi*(static_cast<Wide>(j)+Wide{0.5})*
                              (static_cast<Wide>(k)+Wide{0.5})/
                              static_cast<Wide>(n));
        }
        out[k]=sum;
    }
    return out;
}

template <typename T>
void verify_dct() {
    const Wide tol=std::same_as<T,float>?Wide{0.002}:Wide{2.e-9};
    for(const std::size_t n:{0U,1U,2U,3U,5U,8U,11U,16U,31U,64U,127U})
        for(std::size_t pattern=0;pattern<5;++pattern) {
            const auto x=signal<T>(n,pattern);
            for(int kind=1;kind<=4;++kind) {
                const auto expected=reference_dct<T>(x,kind);
                std::vector<T> direct,via_fft;
                switch(kind) {
                case 1:
                    direct=dct::type1<T>(x);
                    via_fft=dct::type1_bluestein<T>(x);
                    break;
                case 2:
                    direct=dct::type2<T>(x);
                    via_fft=dct::type2_bluestein<T>(x);
                    break;
                case 3:
                    direct=dct::type3<T>(x);
                    via_fft=dct::type3_bluestein<T>(x);
                    break;
                default:
                    direct=dct::type4<T>(x);
                    via_fft=dct::type4_bluestein<T>(x);
                    break;
                }
                assert(direct.size()==n && via_fft.size()==n);
                for(std::size_t k=0;k<n;++k) {
                    check(widen(direct[k]),{expected[k],Wide{}},tol,
                          "dct_direct",n,static_cast<std::size_t>(kind),pattern,k);
                    check(widen(via_fft[k]),{expected[k],Wide{}},tol,
                          "dct_bluestein",n,static_cast<std::size_t>(kind),pattern,k);
                }
            }
        }
}

template <typename Sample>
std::vector<Sample> samples(std::size_t n) {
    using T=signal_processing::detail::scalar_t<Sample>;
    std::vector<Sample> out(n);
    for(std::size_t i=0;i<n;++i) {
        const T re=static_cast<T>(Wide{0.2}*static_cast<Wide>(static_cast<int>(i%9)-4));
        if constexpr (std::same_as<Sample,T>)out[i]=re;
        else out[i]=Sample{re,static_cast<T>(
            Wide{0.15}*static_cast<Wide>(static_cast<int>((3*i)%7)-3))};
    }
    return out;
}

template <typename Sample>
void verify_fir() {
    using T=signal_processing::detail::scalar_t<Sample>;
    const Wide tol=std::same_as<T,float>?Wide{6.e-4}:Wide{5.e-10};
    for(const std::size_t n:{0U,1U,2U,8U,31U,49U})
        for(const std::size_t taps:{0U,1U,2U,3U,7U,9U}) {
            const auto x=samples<Sample>(n);
            const auto h=samples<Sample>(taps);
            std::vector<ComplexWide> expected(n,ComplexWide{});
            for(std::size_t i=0;i<n;++i)
                for(std::size_t j=0;j<=i && j<taps;++j)
                    expected[i]+=widen(h[j])*widen(x[i-j]);
            const auto direct=fir::direct<Sample>(x,h);
            const auto freq=fir::fft<Sample>(x,h);
            assert(direct.size()==n && freq.size()==n);
            for(std::size_t i=0;i<n;++i) {
                check(widen(direct[i]),expected[i],tol,
                      "fir_direct",n,taps,0,i);
                check(widen(freq[i]),expected[i],tol*Wide{4},
                      "fir_fft",n,taps,0,i);
            }
            if(taps==0)continue; // streaming and overlap-save reject empty taps
            const auto overlap=fir::overlap_save<Sample>(x,h,16);
            fir::Filter<Sample> streaming(h);
            const auto streamed=streaming.process(x);
            assert(overlap.size()==n && streamed.size()==n);
            for(std::size_t i=0;i<n;++i) {
                check(widen(overlap[i]),expected[i],tol*Wide{4},
                      "fir_overlap_save",n,taps,0,i);
                check(widen(streamed[i]),expected[i],tol,
                      "fir_streamed",n,taps,0,i);
            }
            streaming.reset();
            for(std::size_t i=0;i<n;++i)
                check(widen(streaming.process(x[i])),expected[i],tol,
                      "fir_reset",n,taps,0,i);
        }
}

template <typename T>
void verify_symmetry_contracts() {
    const T nan=std::numeric_limits<T>::quiet_NaN();
    const T inf=std::numeric_limits<T>::infinity();
    const std::array<T,3> symmetric{T{1},T{2},T{1}};
    const std::array<T,3> antisymmetric{T{1},T{0},T{-1}};
    assert(fir::symmetric_coefficients<T>(symmetric));
    assert(fir::antisymmetric_coefficients<T>(antisymmetric));
    for(const T invalid:{nan,inf,-inf}) {
        // An invalid center tap can otherwise be misclassified as valid.
        const std::array<T,3> bad_center{T{1},invalid,T{1}};
        const std::array<T,3> bad_edge{invalid,T{0},-invalid};
        assert(!fir::symmetric_coefficients<T>(bad_center));
        assert(!fir::antisymmetric_coefficients<T>(bad_center));
        assert(!fir::symmetric_coefficients<T>(bad_edge));
        assert(!fir::antisymmetric_coefficients<T>(bad_edge));
        const std::array<T,2> x{T{1},T{2}};
        bool threw=false;
        try {(void)fir::symmetric_direct<T>(x,bad_center);}
        catch(const std::invalid_argument&) {threw=true;}
        assert(threw);
        threw=false;
        try {(void)fir::antisymmetric_direct<T>(x,bad_edge);}
        catch(const std::invalid_argument&) {threw=true;}
        assert(threw);
    }
    for(const T invalid:{T{-1},nan,inf}) {
        bool threw=false;
        try {(void)fir::symmetric_coefficients<T>(symmetric,invalid);}
        catch(const std::invalid_argument&) {threw=true;}
        assert(threw);
        threw=false;
        try {(void)fir::antisymmetric_coefficients<T>(antisymmetric,invalid);}
        catch(const std::invalid_argument&) {threw=true;}
        assert(threw);
    }
    const std::array<T,5> x{T{1},T{2},T{3},T{4},T{5}};
    const auto ref_sym=fir::direct<T>(x,symmetric);
    const auto ref_asym=fir::direct<T>(x,antisymmetric);
    const auto got_sym=fir::symmetric_direct<T>(x,symmetric);
    const auto got_asym=fir::antisymmetric_direct<T>(x,antisymmetric);
    for(std::size_t i=0;i<x.size();++i) {
        check(widen(got_sym[i]),widen(ref_sym[i]),Wide{1.e-6},
              "linear_phase_symmetric",x.size(),symmetric.size(),0,i);
        check(widen(got_asym[i]),widen(ref_asym[i]),Wide{1.e-6},
              "linear_phase_antisymmetric",x.size(),antisymmetric.size(),0,i);
    }
}
}  // namespace

int main() {
    verify_dct<float>();
    verify_dct<double>();
    verify_fir<float>();
    verify_fir<double>();
    verify_fir<std::complex<float>>();
    verify_fir<std::complex<double>>();
    verify_symmetry_contracts<float>();
    verify_symmetry_contracts<double>();
}
