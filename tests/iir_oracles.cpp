// Independent IIR response oracles with analytically placed conjugate poles.
// Deliberately avoids the library's frequency-response/transfer conversion
// when computing expected outputs.
#ifdef NDEBUG
#error "IIR numerical oracles must have assertions enabled in Release"
#endif

#include "signal_processing/iir.hpp"
#include "signal_processing/iir_design.hpp"

#include <array>
#include <cassert>
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
#include <vector>

namespace {
namespace iir = signal_processing::iir;
namespace design = signal_processing::iir::design;
using L = long double;
using CL = std::complex<L>;

template <iir::Scalar T>
void close_value(T actual, L expected, L tolerance) {
    assert(std::isfinite(actual));
    assert(std::abs(static_cast<L>(actual) - expected) <=
           tolerance * (L{1} + std::abs(expected)));
}

template <iir::Scalar T>
std::vector<T> stimulus(std::size_t count, unsigned pattern, std::uint32_t seed) {
    std::vector<T> input(count);
    std::mt19937 engine(seed);
    std::uniform_real_distribution<T> uniform(T{-1}, T{1});
    for (std::size_t i = 0; i < count; ++i) {
        switch (pattern) {
        case 0: input[i] = i == 0 ? T{1} : T{}; break;  // impulse
        case 1: input[i] = T{1}; break;                 // step
        case 2: input[i] = i % 2 == 0 ? T{1} : T{-1}; break;
        default: input[i] = uniform(engine); break;
        }
    }
    return input;
}

// Textbook difference equation in long double, entirely independent from
// all C++ float/double realization state graphs. Uses the *rounded* public
// coefficient values to distinguish arithmetic error from coefficient error.
template <iir::Scalar T>
std::vector<L> reference_response(
    std::span<const T> input, T b0, T a1, T a2) {
    std::vector<L> output(input.size());
    for (std::size_t n = 0; n < input.size(); ++n) {
        L y = static_cast<L>(b0) * static_cast<L>(input[n]);
        if (n >= 1) y -= static_cast<L>(a1) * output[n - 1];
        if (n >= 2) y -= static_cast<L>(a2) * output[n - 2];
        output[n] = y;
    }
    return output;
}

template <iir::Scalar T>
void verify_one_case(L radius, L angle, unsigned pattern) {
    // Complex-conjugate pole pair p = r exp(±i theta), 0 < r < 1.
    // H(z) = b0 / (1 - 2r cos(theta) z^-1 + r² z^-2);
    // choose b0 so that H(1)=1, without using production design helpers.
    const T r = static_cast<T>(radius);
    const T theta = static_cast<T>(angle);
    const T a1 = static_cast<T>(-2.L * static_cast<L>(r) *
                                std::cos(static_cast<L>(theta)));
    const T a2 = static_cast<T>(static_cast<L>(r) * static_cast<L>(r));
    const T b0 = T{1} + a1 + a2;
    const std::vector<T> numerator{b0};
    const std::vector<T> denominator{T{1},a1,a2};
    const iir::BiquadCoefficients<T> biquad{b0,T{},T{},T{1},a1,a2};

    auto input=stimulus<T>(512,pattern,23711U+pattern*73U+
                                      static_cast<std::uint32_t>(radius*1000.L));
    const auto ref=reference_response<T>(input,b0,a1,a2);
    iir::DirectFormI<T> df1(numerator,denominator);
    iir::DirectFormII<T> df2(numerator,denominator);
    iir::TransposedDirectFormI<T> tdf1(numerator,denominator);
    iir::TransposedDirectFormII<T> tdf2(numerator,denominator);
    iir::BiquadDirectFormI<T> bdf1(biquad);
    iir::BiquadDirectFormII<T> bdf2(biquad);
    iir::BiquadTransposedDirectFormII<T> btdf2(biquad);
    iir::SosCascade<T> cascade({biquad});

    const L tol=std::same_as<T,float> ? 3.e-3L : 2.e-9L;
    for(std::size_t i=0;i<input.size();++i) {
        const T sample=input[i];
        for(const T result:{
            df1.process(sample),df2.process(sample),tdf1.process(sample),
            tdf2.process(sample),bdf1.process(sample),bdf2.process(sample),
            btdf2.process(sample),cascade.process(sample)}) {
            close_value(result,ref[i],tol);
        }
    }
    // Stateful APIs must replay precisely after reset; block processing
    // must also reproduce scalar processing and reject too-small output.
    df1.reset();
    std::vector<T> replay(input.size());
    df1.process(std::span<const T>{input.data(),input.size()},
                std::span<T>{replay.data(),replay.size()});
    for(std::size_t i=0;i<input.size();++i) close_value(replay[i],ref[i],tol);

    bool threw=false;
    try {
        std::array<T,1> insufficient{};
        df2.process(std::span<const T>{input.data(),input.size()},
                    std::span<T>{insufficient.data(),insufficient.size()});
    } catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    df2.reset();
    for(std::size_t i=0;i<input.size();++i)
        close_value(df2.process(input[i]),ref[i],tol);

    // A separate z-plane mathematical oracle checks the public frequency
    // response over DC, Nyquist, and several interior frequencies.
    const T real=static_cast<T>(static_cast<L>(r)*std::cos(static_cast<L>(theta)));
    const T imag=static_cast<T>(static_cast<L>(r)*std::sin(static_cast<L>(theta)));
    const design::Zpk<T> zpk{{T{},T{}},
                             {{real,imag},{real,-imag}},b0};
    assert(design::stable(zpk));
    for(const T frequency:{T{},static_cast<T>(0.075),static_cast<T>(0.3),
                           static_cast<T>(0.7),T{1}}) {
        const L omega=std::numbers::pi_v<L>*static_cast<L>(frequency);
        const CL zinv{std::cos(omega),-std::sin(omega)};
        const CL den=CL{1,0}+static_cast<L>(a1)*zinv+
                     static_cast<L>(a2)*zinv*zinv;
        const CL expected=static_cast<L>(b0)/den;
        const auto actual=design::frequency_response(zpk,frequency);
        const CL observed{static_cast<L>(actual.real()),
                          static_cast<L>(actual.imag())};
        assert(std::isfinite(actual.real())&&std::isfinite(actual.imag()));
        const L limit=std::same_as<T,float> ? 4.e-3L : 3.e-10L;
        assert(std::abs(observed-expected)<=limit*(1.L+std::abs(expected)));
    }
}

// Closed-form independent oracle: H(z)=(1-r)/(1-r z^-1) has impulse
// response h[n]=(1-r)r^n. Includes negative poles and non-unit a0.
template <iir::Scalar T>
void verify_first_order_impulse() {
    const L tolerance=std::same_as<T,float> ? 1.e-4L : 3.e-12L;
    for(const L pole:{-0.93L,-0.25L,0.17L,0.87L,0.985L}) {
        const T r=static_cast<T>(pole);
        const T b0=T{1}-r;
        const iir::FirstOrderCoefficients<T> c{
            T{2}*b0,T{},T{2},-T{2}*r};
        iir::FirstOrderTransposedDirectFormII<T> first(c);
        iir::DirectFormI<T> general({b0},{T{1},-r});
        L expected=static_cast<L>(b0);
        for(std::size_t n=0;n<384;++n) {
            const T input=n==0 ? T{1} : T{};
            close_value(first.process(input),expected,tolerance);
            close_value(general.process(input),expected,tolerance);
            expected *= static_cast<L>(r);
        }
        first.reset();
        close_value(first.process(T{1}),static_cast<L>(b0),tolerance);
    }
}

template <iir::Scalar T>
void invalid_domains() {
    bool threw=false;
    try { (void)iir::DirectFormI<T>({T{1}},{T{0},T{1}}); }
    catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    threw=false;
    try { (void)iir::DirectFormII<T>({}, {T{1}}); }
    catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    threw=false;
    try { (void)iir::BiquadTransposedDirectFormII<T>(
        {T{1},T{},T{},T{},T{},T{}}); }
    catch(const std::invalid_argument&) {threw=true;}
    assert(threw);
    const design::Zpk<T> unstable{{},{{static_cast<T>(1.01),T{}}},T{1}};
    assert(!design::stable(unstable));
    const design::Zpk<T> unit{{},{{T{1},T{}}},T{1}};
    assert(!design::stable(unit));
    const design::Zpk<T> stable{{},{{static_cast<T>(0.95),T{}}},T{1}};
    for(const T invalid:{T{-1},static_cast<T>(1.01),
                         std::numeric_limits<T>::quiet_NaN()}) {
        threw=false;
        try { (void)design::frequency_response(stable,invalid); }
        catch(const std::invalid_argument&) {threw=true;}
        assert(threw);
    }
}

template <iir::Scalar T>
void run() {
    invalid_domains<T>();
    verify_first_order_impulse<T>();
    for(const L radius:{0.05L,0.43L,0.82L,0.985L})
        for(const L angle:{0.13L,0.47L,1.12L,2.37L})
            for(unsigned pattern=0;pattern<4;++pattern)
                verify_one_case<T>(radius,angle,pattern);
}
}  // namespace

int main() {
    run<float>();
    run<double>();
}
