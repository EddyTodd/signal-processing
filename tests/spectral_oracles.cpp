// Independent, long-double textbook STFT / density oracles.
// The expected spectra do not call production FFT or spectral algorithms.
#ifdef NDEBUG
#error "Spectral numerical oracles must run with assertions enabled"
#endif

#include "signal_processing/stft.hpp"
#include "signal_processing/spectral.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <span>
#include <stdexcept>
#include <vector>

namespace {
namespace stft = signal_processing::stft;
namespace spectral = signal_processing::spectral;
using Wide = long double;
using CW = std::complex<Wide>;

template <typename Sample>
using Scalar = signal_processing::detail::scalar_t<Sample>;

template <typename Sample>
CW widen(Sample sample) {
    if constexpr (std::same_as<Sample,float> || std::same_as<Sample,double>) {
        return {static_cast<Wide>(sample), Wide{}};
    } else {
        return {static_cast<Wide>(sample.real()),
                static_cast<Wide>(sample.imag())};
    }
}

void check_close(CW actual, CW expected, Wide tolerance,
                 const char* name, std::size_t n, std::size_t frame,
                 std::size_t hop, std::size_t index) {
    const Wide error=std::abs(actual-expected);
    const Wide bound=tolerance*(Wide{1}+std::abs(expected));
    if(!std::isfinite(actual.real())||!std::isfinite(actual.imag())||
       !(error<=bound)) {
        std::cerr << "SPECTRAL_ORACLE_FAILURE family=" << name
                  << " n=" << n << " frame=" << frame << " hop=" << hop
                  << " index=" << index << " actual=" << actual
                  << " expected=" << expected << " error=" << error
                  << " bound=" << bound << '\n';
        std::abort();
    }
}

template <typename Sample>
std::vector<Sample> signal(std::size_t size, Wide phase=Wide{}) {
    using T=Scalar<Sample>;
    std::vector<Sample> result(size);
    for(std::size_t i=0;i<size;++i) {
        const Wide x=static_cast<Wide>(i);
        const T re=static_cast<T>(0.37L*std::sin(0.43L*x+phase) +
                                  0.11L*std::cos(0.67L*x) +
                                  0.021L*static_cast<Wide>(i%3));
        if constexpr (std::same_as<Sample,T>) {
            result[i]=re;
        } else {
            const T im=static_cast<T>(0.29L*std::cos(0.31L*x+phase) -
                                      0.13L*std::sin(0.83L*x));
            result[i]=Sample{re,im};
        }
    }
    return result;
}

template <typename T>
std::vector<T> window(std::size_t size) {
    std::vector<T> w(size);
    for(std::size_t i=0;i<size;++i)
        w[i]=static_cast<T>(0.5L+0.05L*static_cast<Wide>(i%5));
    return w;
}

// Independent enumeration of signal-bearing frame starts. With a hop
// longer than the window, gaps are allowed; no empty frame past n is emitted.
std::vector<std::size_t> starts(std::size_t n, std::size_t length,
                                std::size_t hop, bool pad) {
    std::vector<std::size_t> result;
    if(n==0)return result;
    if(!pad) {
        if(n<length)return result;
        for(std::size_t offset=0;offset<=n-length;offset+=hop)
            result.push_back(offset);
        return result;
    }
    for(std::size_t offset=0;offset<n;offset+=hop) {
        result.push_back(offset);
        if(length>=n-offset)break;
    }
    return result;
}

// Direct O(L^2) DFT in long double, with explicit zero padding.
template <typename Sample>
CW reference_bin(const std::vector<Sample>& x,
                 const std::vector<Scalar<Sample>>& w,
                 std::size_t offset, std::size_t k) {
    CW sum{};
    for(std::size_t j=0;j<w.size();++j) {
        if(offset+j>=x.size())continue;
        const Wide angle=-Wide{2}*std::numbers::pi_v<Wide>*
                          static_cast<Wide>(j)*static_cast<Wide>(k)/
                          static_cast<Wide>(w.size());
        sum+=widen(x[offset+j])*static_cast<Wide>(w[j])*
             CW{std::cos(angle),std::sin(angle)};
    }
    return sum;
}

template <typename Sample>
void verify_stft(std::size_t n, std::size_t length,
                 std::size_t hop, bool pad) {
    using T=Scalar<Sample>;
    const auto x=signal<Sample>(n);
    const auto w=window<T>(length);
    const auto positions=starts(n,length,hop,pad);
    const auto actual=stft::bluestein<Sample>(x,w,hop,pad);
    assert(actual.original_size==n);
    assert(actual.frame_size==length && actual.hop_size==hop);
    assert(actual.frames.size()==positions.size());
    const Wide tol=std::same_as<T,float> ? 8.e-4L : 5.e-10L;
    for(std::size_t m=0;m<positions.size();++m) {
        assert(actual.frames[m].size()==length);
        for(std::size_t k=0;k<length;++k)
            check_close(widen(actual.frames[m][k]),
                        reference_bin(x,w,positions[m],k),
                        tol,"stft_dft",n,length,hop,m*length+k);
    }
    // The inverse must reconstruct each sample with nonzero overlap
    // denominator. Uncovered gaps and truncated tails explicitly yield zero.
    const auto reconstructed=stft::inverse_bluestein<Sample>(actual,w);
    assert(reconstructed.size()==n);
    for(std::size_t i=0;i<n;++i) {
        bool covered=false;
        for(const std::size_t offset:positions)
            if(i>=offset && i-offset<length)covered=true;
        const CW expected=covered?widen(x[i]):CW{};
        check_close(widen(reconstructed[i]),expected,tol*Wide{4},
                    "istft_overlap",n,length,hop,i);
    }
    const auto power=stft::power_spectrogram<Sample>(actual,false);
    assert(power.size()==positions.size());
    for(std::size_t m=0;m<power.size();++m) {
        assert(power[m].size()==length);
        for(std::size_t k=0;k<length;++k) {
            const Wide expected=std::norm(reference_bin(x,w,positions[m],k));
            check_close(CW{static_cast<Wide>(power[m][k]),0},
                        CW{expected,0},tol*Wide{8},
                        "stft_power",n,length,hop,m*length+k);
        }
    }
    if constexpr (std::same_as<Sample,T>) {
        const auto one=stft::power_spectrogram<Sample>(actual,true);
        for(const auto& frame:one)assert(frame.size()==length/2+1);
    } else {
        bool threw=false;
        try { (void)stft::power_spectrogram<Sample>(actual,true); }
        catch(const std::invalid_argument&) {threw=true;}
        assert(threw);
    }
}

template <typename Sample>
void verify_density(std::size_t length, std::size_t hop,
                    spectral::Sides sides) {
    using T=Scalar<Sample>;
    constexpr std::size_t n=19;
    const auto x=signal<Sample>(n);
    const auto y=signal<Sample>(n,Wide{0.71});
    const auto w=window<T>(length);
    const T rate=T{37};
    Wide energy{};
    for(const T value:w)energy+=static_cast<Wide>(value)*static_cast<Wide>(value);
    const Wide scale=Wide{1}/(static_cast<Wide>(rate)*energy);
    const std::size_t count=sides==spectral::Sides::one_sided?length/2+1:length;
    const auto positions=starts(n,length,hop,false);
    const Wide tol=std::same_as<T,float> ? 1.e-3L : 1.e-9L;
    const auto psd=spectral::periodogram_bluestein<Sample>(
        std::span<const Sample>(x.data(),length),w,rate,sides);
    assert(psd.size()==count);
    for(std::size_t k=0;k<count;++k) {
        Wide value=std::norm(reference_bin(x,w,0,k))*scale;
        if(sides==spectral::Sides::one_sided && k!=0 &&
           !(length%2==0 && k==length/2))value*=Wide{2};
        check_close(CW{static_cast<Wide>(psd[k]),0},CW{value,0},tol,
                    "periodogram",n,length,hop,k);
    }
    const auto welch=spectral::welch_bluestein<Sample>(x,w,hop,rate,sides);
    const auto csd=spectral::cross_spectral_density_bluestein<Sample>(
        x,y,w,hop,rate,sides);
    const auto gram=spectral::spectrogram_bluestein<Sample>(
        x,w,hop,rate,sides,false);
    assert(welch.size()==count && csd.size()==count);
    assert(gram.size()==positions.size());
    for(std::size_t k=0;k<count;++k) {
        Wide sum{};
        CW cross{};
        for(std::size_t m=0;m<positions.size();++m) {
            const auto left=reference_bin(x,w,positions[m],k);
            const auto right=reference_bin(y,w,positions[m],k);
            const Wide multiplier=(sides==spectral::Sides::one_sided &&
                                   k!=0 && !(length%2==0 && k==length/2))
                                  ? Wide{2} : Wide{1};
            const Wide expected=std::norm(left)*scale*multiplier;
            assert(gram[m].size()==count);
            check_close(CW{static_cast<Wide>(gram[m][k]),0},
                        CW{expected,0},tol,"spectrogram",n,length,hop,
                        m*length+k);
            sum+=expected;
            cross+=left*std::conj(right)*scale*multiplier;
        }
        const Wide inv=Wide{1}/static_cast<Wide>(positions.size());
        check_close(CW{static_cast<Wide>(welch[k]),0},
                    CW{sum*inv,0},tol,"welch",n,length,hop,k);
        check_close(widen(csd[k]),cross*inv,tol,
                    "cross_density",n,length,hop,k);
    }
}

template <typename Sample>
void verify_invalid_contracts() {
    using T=Scalar<Sample>;
    const auto x=signal<Sample>(8);
    const auto w=window<T>(5);
    bool threw=false;
    try { (void)stft::bluestein<Sample>(x,std::span<const T>{},1,true); }
    catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    threw=false;
    try { (void)stft::bluestein<Sample>(x,w,0,true); }
    catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    auto transform=stft::bluestein<Sample>(x,w,2,true);
    transform.frames[0].pop_back();
    threw=false;
    try { (void)stft::inverse_bluestein<Sample>(transform,w); }
    catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    threw=false;
    try { (void)spectral::periodogram_bluestein<Sample>(
        std::span<const Sample>(x.data(),w.size()),w,T{},spectral::Sides::two_sided); }
    catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    const std::vector<T> zeros(w.size(),T{});
    threw=false;
    try { (void)spectral::periodogram_bluestein<Sample>(
        std::span<const Sample>(x.data(),w.size()),zeros,T{37},
        spectral::Sides::two_sided); }
    catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    if constexpr (!std::same_as<Sample,T>) {
        threw=false;
        try { (void)spectral::periodogram_bluestein<Sample>(
            std::span<const Sample>(x.data(),w.size()),w,T{37},
            spectral::Sides::one_sided); }
        catch(const std::invalid_argument&) {threw=true;}
        assert(threw);
    }
}

template <typename Sample>
void run() {
    for(const std::size_t n:{0U,1U,2U,5U,9U,11U,17U})
        for(const std::size_t length:{3U,4U,8U})
            for(const std::size_t hop:{1U,2U,4U,9U})
                for(const bool pad:{false,true})
                    verify_stft<Sample>(n,length,hop,pad);
    for(const std::size_t length:{4U,5U,8U})
        for(const std::size_t hop:{2U,3U,8U}) {
            verify_density<Sample>(length,hop,spectral::Sides::two_sided);
            if constexpr (std::same_as<Sample,Scalar<Sample>>)
                verify_density<Sample>(length,hop,spectral::Sides::one_sided);
        }
    verify_invalid_contracts<Sample>();
}
}  // namespace

int main() {
    run<float>();
    run<double>();
    run<std::complex<float>>();
    run<std::complex<double>>();
}
