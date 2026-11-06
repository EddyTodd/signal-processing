// Strict spectral parameter-domain regression, including early returns.
#ifdef NDEBUG
#error "Spectral domain tests require Release assertions"
#endif
#include "signal_processing/spectral.hpp"
#include <cassert>
#include <complex>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace {
namespace spectral=signal_processing::spectral;

template<class F>
void must_reject(F&& f) {
    bool rejected=false;
    try { (void)f(); }
    catch(const std::invalid_argument&) { rejected=true; }
    assert(rejected);
}

template<class Sample>
void verify() {
    using T=signal_processing::detail::scalar_t<Sample>;
    const std::vector<Sample> signal(8,Sample{1});
    const std::vector<T> window(5,T{1});
    const std::span<const Sample> segment(signal.data(),window.size());
    constexpr T rate=T{37};
    const auto reject_every_method=[&](std::span<const T> win,T frequency,
                                       spectral::Sides side) {
        must_reject([&]{return spectral::periodogram_bluestein<Sample>(
            segment,win,frequency,side);});
        must_reject([&]{return spectral::welch_bluestein<Sample>(
            signal,win,2,frequency,side);});
        must_reject([&]{return spectral::cross_spectral_density_bluestein<Sample>(
            signal,signal,win,2,frequency,side);});
        must_reject([&]{return spectral::spectrogram_bluestein<Sample>(
            signal,win,2,frequency,side);});
    };
    const T nan=std::numeric_limits<T>::quiet_NaN();
    const T inf=std::numeric_limits<T>::infinity();
    for(const T invalid:{nan,inf,-inf,T{},T{-1}})
        reject_every_method(window,invalid,spectral::Sides::two_sided);
    reject_every_method(window,rate,static_cast<spectral::Sides>(99));
    const std::vector<T> zeros(5,T{});
    reject_every_method(zeros,rate,spectral::Sides::two_sided);
    for(const T invalid:{nan,inf,-inf}) {
        auto bad=window;
        bad[2]=invalid;
        reject_every_method(bad,rate,spectral::Sides::two_sided);
    }
    reject_every_method(
        std::vector<T>(5,std::numeric_limits<T>::max()/T{2}),
        rate,spectral::Sides::two_sided);
    reject_every_method(window,std::numeric_limits<T>::max(),
                        spectral::Sides::two_sided);
    reject_every_method(window,std::numeric_limits<T>::denorm_min(),
                        spectral::Sides::two_sided);

    const std::vector<Sample> empty;
    assert(spectral::periodogram_bluestein<Sample>(
        empty,std::span<const T>{},rate,spectral::Sides::two_sided).empty());
    must_reject([&]{return spectral::periodogram_bluestein<Sample>(
        empty,window,rate,spectral::Sides::two_sided);});
    // A return with no full frames must not mask invalid window energy.
    const auto short_signal=segment.subspan(0,2);
    must_reject([&]{return spectral::welch_bluestein<Sample>(
        short_signal,zeros,2,rate,spectral::Sides::two_sided);});
    must_reject([&]{return spectral::cross_spectral_density_bluestein<Sample>(
        short_signal,short_signal,zeros,2,rate,spectral::Sides::two_sided);});
    // The same undersized input with a finite window remains an empty result.
    assert(spectral::welch_bluestein<Sample>(
        short_signal,window,2,rate,spectral::Sides::two_sided).empty());
    assert(spectral::cross_spectral_density_bluestein<Sample>(
        short_signal,short_signal,window,2,rate,
        spectral::Sides::two_sided).empty());
}
}
int main() {
    verify<float>();
    verify<double>();
    verify<std::complex<float>>();
    verify<std::complex<double>>();
}
