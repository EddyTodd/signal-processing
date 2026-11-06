// Independent long-double mathematical linear/circular convolution and
// cross-correlation oracle, including partitioned and streaming state graphs.
#ifdef NDEBUG
#error "Convolution oracles require enabled assertions in Release builds"
#endif

#include "signal_processing/convolution.hpp"
#include "signal_processing/correlation.hpp"

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
#include <span>
#include <stdexcept>
#include <vector>

namespace {
namespace conv=signal_processing::convolution;
namespace corr=signal_processing::correlation;
using Wide=long double;
using CW=std::complex<Wide>;

template <typename Sample>
using Scalar=signal_processing::detail::scalar_t<Sample>;

template <typename Sample>
CW to_wide(Sample x) {
    if constexpr (std::same_as<Sample,Scalar<Sample>>)
        return {static_cast<Wide>(x),Wide{}};
    else
        return {static_cast<Wide>(x.real()),static_cast<Wide>(x.imag())};
}

void require_close(CW actual,CW expected,Wide tol,
                   const char* name,std::size_t n,std::size_t m,
                   std::size_t shape,std::size_t index) {
    const Wide err=std::abs(actual-expected);
    const Wide bound=tol*(Wide{1}+std::abs(expected));
    if(!std::isfinite(actual.real())||!std::isfinite(actual.imag())||
       !(err<=bound)) {
        std::cerr<<"CONV_ORACLE_FAILURE method="<<name
                 <<" n="<<n<<" m="<<m<<" shape="<<shape
                 <<" index="<<index<<" actual="<<actual
                 <<" expected="<<expected<<" abs_error="<<err
                 <<" bound="<<bound<<'\n';
        std::abort();
    }
}

template <typename Sample>
std::vector<Sample> signal(std::size_t n,std::size_t pattern) {
    using T=Scalar<Sample>;
    std::vector<Sample> result(n);
    for(std::size_t i=0;i<n;++i) {
        const Wide re=pattern==0
            ?Wide{0.2}*static_cast<Wide>(static_cast<int>((11*i+3)%13)-6)
            :pattern==1
                ?(i==0?Wide{1}:Wide{})
                :(i%3==0?Wide{0}:Wide{0.035}*static_cast<Wide>(
                    static_cast<int>((7*i+5)%11)-5));
        const T r=static_cast<T>(re);
        if constexpr (std::same_as<Sample,T>)result[i]=r;
        else result[i]=Sample{r,static_cast<T>(
            pattern==1?Wide{}:Wide{0.11}*static_cast<Wide>(
                static_cast<int>((3*i+1)%7)-3))};
    }
    return result;
}

template <typename Sample>
std::vector<CW> linear_reference(const std::vector<Sample>& a,
                                 const std::vector<Sample>& b) {
    if(a.empty()||b.empty())return {};
    std::vector<CW> out(a.size()+b.size()-1,CW{});
    for(std::size_t i=0;i<a.size();++i)
        for(std::size_t j=0;j<b.size();++j)
            out[i+j]+=to_wide(a[i])*to_wide(b[j]);
    return out;
}

template <typename Sample>
std::vector<CW> cross_reference(const std::vector<Sample>& a,
                                const std::vector<Sample>& b,
                                bool normalized) {
    if(a.empty()||b.empty())return {};
    std::vector<CW> out(a.size()+b.size()-1,CW{});
    for(std::size_t k=0;k<out.size();++k) {
        const auto lag=static_cast<std::ptrdiff_t>(k)-
                       static_cast<std::ptrdiff_t>(b.size()-1);
        Wide ea{},eb{};
        for(std::size_t i=0;i<a.size();++i) {
            const auto j=static_cast<std::ptrdiff_t>(i)-lag;
            if(j<0||j>=static_cast<std::ptrdiff_t>(b.size()))continue;
            const CW x=to_wide(a[i]),y=to_wide(b[static_cast<std::size_t>(j)]);
            out[k]+=x*std::conj(y);
            ea+=std::norm(x);
            eb+=std::norm(y);
        }
        if(normalized) {
            const Wide denom=std::sqrt(ea)*std::sqrt(eb);
            out[k]=denom>Wide{}?out[k]/denom:CW{};
        }
    }
    return out;
}

template <typename Sample>
void compare(const std::vector<Sample>& actual,const std::vector<CW>& reference,
             Wide tol,const char* name,std::size_t n,std::size_t m,
             std::size_t shape=0) {
    assert(actual.size()==reference.size());
    for(std::size_t i=0;i<reference.size();++i)
        require_close(to_wide(actual[i]),reference[i],tol,name,n,m,shape,i);
}

template <typename Sample>
void linear_case(std::size_t n,std::size_t m,std::size_t pattern) {
    const auto a=signal<Sample>(n,pattern),b=signal<Sample>(m,2-pattern);
    const auto expected=linear_reference(a,b);
    const Wide tol=std::same_as<Scalar<Sample>,float>?Wide{0.0004}:Wide{3.e-10};
    compare(conv::direct<Sample>(a,b),expected,tol,"direct",n,m,pattern);
    compare(conv::fft<Sample>(a,b),expected,tol*Wide{5},"fft",n,m,pattern);
    compare(corr::cross_direct<Sample>(a,b),
            cross_reference(a,b,false),tol,"cross_direct",n,m,pattern);
    compare(corr::cross_fft<Sample>(a,b),
            cross_reference(a,b,false),tol*Wide{5},"cross_fft",n,m,pattern);
    compare(corr::normalized_direct<Sample>(a,b),
            cross_reference(a,b,true),tol*Wide{5},
            "normalized_direct",n,m,pattern);
    compare(corr::normalized_fft<Sample>(a,b),
            cross_reference(a,b,true),tol*Wide{7},
            "normalized_fft",n,m,pattern);
    if(n==0||m==0)return;
    for(const std::size_t block:{1U,3U,8U}) {
        compare(conv::overlap_add<Sample>(a,b,block),expected,
                tol*Wide{7},"overlap_add",n,m,block);
        compare(conv::partitioned<Sample>(a,b,block),expected,
                tol*Wide{10},"partitioned",n,m,block);
    }
    for(const std::size_t transform:{8U,16U,32U}) {
        if(m>transform)continue;
        compare(conv::overlap_save<Sample>(a,b,transform),expected,
                tol*Wide{10},"overlap_save",n,m,transform);
    }
    conv::StreamingDirect<Sample> direct_stream(b);
    auto streamed=direct_stream.process(a);
    const auto tail=direct_stream.flush();
    streamed.insert(streamed.end(),tail.begin(),tail.end());
    compare(streamed,expected,tol,"streaming_flush",n,m,pattern);
    direct_stream.reset();
    auto replay=direct_stream.process(a);
    auto replay_tail=direct_stream.flush();
    replay.insert(replay.end(),replay_tail.begin(),replay_tail.end());
    compare(replay,expected,tol,"streaming_reset",n,m,pattern);

    for(const std::size_t block:{1U,3U,8U}) {
        conv::UniformPartitionedConvolver<Sample> stream(b,block);
        std::vector<Sample> in(block,Sample{}),out(block,Sample{}),joined;
        for(std::size_t offset=0;joined.size()<expected.size();offset+=block) {
            std::fill(in.begin(),in.end(),Sample{});
            if(offset<a.size()) {
                const std::size_t count=std::min(block,a.size()-offset);
                std::copy_n(a.begin()+static_cast<std::ptrdiff_t>(offset),
                            count,in.begin());
            }
            stream.process_block(in,out);
            const std::size_t keep=std::min(block,expected.size()-joined.size());
            joined.insert(joined.end(),out.begin(),
                          out.begin()+static_cast<std::ptrdiff_t>(keep));
        }
        compare(joined,expected,tol*Wide{10},
                "partitioned_blocks",n,m,block);
        stream.reset();
        std::fill(in.begin(),in.end(),Sample{});
        std::fill(out.begin(),out.end(),Sample{});
        std::copy_n(a.begin(),std::min(block,a.size()),in.begin());
        stream.process_block(in,out);
        for(std::size_t i=0;i<std::min(block,expected.size());++i)
            require_close(to_wide(out[i]),expected[i],tol*Wide{10},
                         "partitioned_reset",n,m,block,i);
    }
}

template <typename Sample>
void verify_sparse_end_energy() {
    using T=Scalar<Sample>;
    for(const Wide weak:{Wide{0.0001},Wide{0.000000001}}) {
        std::vector<Sample> x(23,Sample{});
        if constexpr (std::same_as<Sample,T>) {
            x.front()=T{1};
            x.back()=static_cast<T>(weak);
        } else {
            x.front()=Sample{T{1},T{}};
            x.back()=Sample{static_cast<T>(weak),T{}};
        }
        const auto expected=cross_reference(x,x,true);
        const auto normal_direct=corr::normalized_direct<Sample>(x,x);
        const auto normal_fft=corr::normalized_fft<Sample>(x,x);
        compare(normal_direct,expected,Wide{0.0005},
                "sparse_tail_normalized_direct",x.size(),x.size());
        compare(normal_fft,expected,Wide{0.0005},
                "sparse_tail_normalized_fft",x.size(),x.size());
        // Both single-sample overlaps must normalize to positive one.
        require_close(to_wide(normal_direct.front()),CW{1,0},Wide{0.0005},
                     "weak_lag",x.size(),x.size(),0,0);
    }
}

template <typename Sample>
void invalid_domains() {
    const std::vector<Sample> x{Sample{1},Sample{2},Sample{3}};
    bool threw=false;
    try{(void)conv::overlap_add<Sample>(x,x,0);}
    catch(const std::invalid_argument&){threw=true;}
    assert(threw);
    threw=false;
    try{(void)conv::OverlapSavePlan<Sample>(x,3);}
    catch(const std::invalid_argument&){threw=true;}
    assert(threw);
    threw=false;
    try{(void)conv::UniformPartitionedConvolver<Sample>(x,0);}
    catch(const std::invalid_argument&){threw=true;}
    assert(threw);
    threw=false;
    try{(void)conv::StreamingDirect<Sample>(std::span<const Sample>{});}
    catch(const std::invalid_argument&){threw=true;}
    assert(threw);
    conv::UniformPartitionedConvolver<Sample> convolver(x,2);
    std::array<Sample,1> bad{},output{};
    threw=false;
    try{convolver.process_block(bad,output);}
    catch(const std::invalid_argument&){threw=true;}
    assert(threw);
}

template <typename Sample>
void run() {
    for(const std::size_t n:{0U,1U,2U,5U,11U,29U})
        for(const std::size_t m:{0U,1U,2U,3U,7U,15U,21U})
            for(const std::size_t pattern:{0U,1U,2U})
                linear_case<Sample>(n,m,pattern);
    verify_sparse_end_energy<Sample>();
    invalid_domains<Sample>();
}
} // namespace

int main(){
    run<float>();
    run<double>();
    run<std::complex<float>>();
    run<std::complex<double>>();
}
