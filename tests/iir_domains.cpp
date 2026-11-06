// Exceptional coefficient domains and independent fourth-order SOS recurrence.
// Keeps all Release assertions active; doesn't relax existing IIR tolerances.
#ifdef NDEBUG
#error "IIR domain and SOS oracles require active Release assertions"
#endif

#include "signal_processing/iir.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace {
namespace iir=signal_processing::iir;
using Wide=long double;

template <typename F>
void rejects(F&& attempt) {
    bool threw=false;
    try { attempt(); }
    catch(const std::invalid_argument&) { threw=true; }
    assert(threw);
}

template <iir::Scalar T>
void invalid_coefficients() {
    const T nan=std::numeric_limits<T>::quiet_NaN();
    const T inf=std::numeric_limits<T>::infinity();
    const T tiny=std::numeric_limits<T>::min();
    const T huge=std::numeric_limits<T>::max();

    rejects([&]{ (void)iir::normalize<T>({T{1}},{nan}); });
    rejects([&]{ (void)iir::normalize<T>({T{1}},{inf}); });
    rejects([&]{ (void)iir::normalize<T>({nan},{T{1}}); });
    rejects([&]{ (void)iir::normalize<T>({T{1},inf},{T{1}}); });
    rejects([&]{ (void)iir::normalize<T>({T{1}},{T{1},-inf}); });
    rejects([&]{ (void)iir::normalize<T>({huge},{tiny}); });
    rejects([&]{ (void)iir::normalize<T>({T{1}},{tiny,huge}); });

    for(const T bad:{nan,inf,-inf}) {
        for(const std::vector<T> numerator : {
                std::vector<T>{bad},std::vector<T>{T{1},bad}}) {
            const std::vector<T> denominator{T{1},T{-0.4}};
            rejects([&]{ (void)iir::DirectFormI<T>(numerator,denominator); });
            rejects([&]{ (void)iir::DirectFormII<T>(numerator,denominator); });
            rejects([&]{ (void)iir::TransposedDirectFormI<T>(numerator,denominator); });
            rejects([&]{ (void)iir::TransposedDirectFormII<T>(numerator,denominator); });
        }
        const iir::FirstOrderCoefficients<T> first{T{1},bad,T{1},T{-0.5}};
        rejects([&]{ (void)iir::FirstOrderTransposedDirectFormII<T>(first); });
        const iir::FirstOrderCoefficients<T> first_a{T{1},T{},bad,T{-0.5}};
        rejects([&]{ (void)iir::FirstOrderTransposedDirectFormII<T>(first_a); });
        const std::array<iir::BiquadCoefficients<T>,3> biquads{{
            {bad,T{},T{},T{1},T{-0.5},T{0.125}},
            {T{1},T{},T{},bad,T{-0.5},T{0.125}},
            {T{1},T{},T{},T{1},T{-0.5},bad}
        }};
        for(const auto c:biquads) {
            rejects([&]{ (void)iir::BiquadDirectFormI<T>(c); });
            rejects([&]{ (void)iir::BiquadDirectFormII<T>(c); });
            rejects([&]{ (void)iir::BiquadTransposedDirectFormII<T>(c); });
            rejects([&]{ (void)iir::SosCascade<T>({c}); });
        }
        rejects([&]{ (void)iir::SosCascade<T>({},bad); });
    }

    const iir::FirstOrderCoefficients<T> first_overflow{
        huge,T{0},tiny,T{0}};
    rejects([&]{ (void)iir::FirstOrderTransposedDirectFormII<T>(first_overflow); });
    const iir::BiquadCoefficients<T> biquad_overflow{
        huge,T{},T{},tiny,T{},T{}};
    rejects([&]{ (void)iir::BiquadDirectFormI<T>(biquad_overflow); });
    rejects([&]{ (void)iir::BiquadDirectFormII<T>(biquad_overflow); });
    rejects([&]{ (void)iir::BiquadTransposedDirectFormII<T>(biquad_overflow); });
    rejects([&]{ (void)iir::SosCascade<T>({biquad_overflow}); });

    // Subnormal and negative finite a0 are legitimate when their
    // normalized ratios remain finite; do not reject by magnitude alone.
    const auto normalized=iir::normalize<T>({tiny,T{0}},{tiny,-tiny});
    assert(normalized.numerator[0]==T{1});
    assert(normalized.denominator[0]==T{1});
    assert(normalized.denominator[1]==T{-1});
    const auto negative=iir::normalize<T>({T{-2},T{1}},{T{-2},T{1}});
    assert(negative.numerator[0]==T{1});
    assert(negative.denominator[0]==T{1});
    assert(negative.numerator[1]==T{-0.5});
}

template <iir::Scalar T>
void finite_fourth_order() {
    const Wide tolerance=std::same_as<T,float>?Wide{0.003}:Wide{8.e-10};
    const std::array<Wide,2> radii{Wide{0.76},Wide{0.985}};
    const std::array<Wide,2> angles{Wide{0.44},Wide{1.13}};
    std::array<iir::BiquadCoefficients<T>,2> sections{};
    for(std::size_t i=0;i<2;++i) {
        const T r=static_cast<T>(radii[i]);
        const T theta=static_cast<T>(angles[i]);
        const T a1=static_cast<T>(-Wide{2}*static_cast<Wide>(r)*
                                  std::cos(static_cast<Wide>(theta)));
        const T a2=static_cast<T>(static_cast<Wide>(r)*static_cast<Wide>(r));
        const T b0=T{1}+a1+a2;
        sections[i]={b0,T{},T{},T{1},a1,a2};
    }
    // Closed recurrence of two independent sections in long double:
    // y_k[n]=b0_k * y_(k-1)[n] - a1_k*y_k[n-1] - a2_k*y_k[n-2].
    // The source coefficients are rounded to public T before widening.
    for(std::size_t pattern=0;pattern<3;++pattern) {
        iir::SosCascade<T> cascade({sections[0],sections[1]});
        iir::BiquadTransposedDirectFormII<T> first(sections[0]),second(sections[1]);
        std::array<Wide,2> y1{},y2{};
        for(std::size_t n=0;n<512;++n) {
            const T stimulus=pattern==0?(n==0?T{1}:T{}):
                pattern==1?T{1}:static_cast<T>(
                    Wide{0.41}*std::sin(Wide{0.29}*static_cast<Wide>(n))+
                    Wide{0.17}*std::cos(Wide{0.71}*static_cast<Wide>(n)));
            const Wide in=static_cast<Wide>(stimulus);
            Wide x=in;
            for(std::size_t stage=0;stage<2;++stage) {
                const auto c=sections[stage];
                const Wide out=static_cast<Wide>(c.b0)*x-
                               static_cast<Wide>(c.a1)*y1[stage]-
                               static_cast<Wide>(c.a2)*y2[stage];
                y2[stage]=y1[stage];
                y1[stage]=out;
                x=out;
            }
            const T actual=cascade.process(stimulus);
            const T separately=second.process(first.process(stimulus));
            const Wide error=std::abs(static_cast<Wide>(actual)-x);
            const Wide bound=tolerance*(Wide{1}+std::abs(x));
            if(!std::isfinite(actual)||!(error<=bound)) {
                std::cerr<<"IIR_DOMAIN_ORACLE_FAILURE precision="
                         <<(std::same_as<T,float>?"float32":"float64")
                         <<" pattern="<<pattern<<" sample="<<n
                         <<" actual="<<actual<<" expected="<<x
                         <<" error="<<error<<" bound="<<bound<<'\n';
                std::abort();
            }
            assert(std::abs(static_cast<Wide>(actual)-static_cast<Wide>(separately))
                   <= tolerance*(Wide{1}+std::abs(x)));
        }
        cascade.reset();
        assert(std::abs(cascade.process(T{1})-sections[0].b0*sections[1].b0)
               <= static_cast<T>(tolerance));
    }
}

template <iir::Scalar T>
void run() {
    invalid_coefficients<T>();
    finite_fourth_order<T>();
}
}  // namespace

int main() {
    run<float>();
    run<double>();
}
