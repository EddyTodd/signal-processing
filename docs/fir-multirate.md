# FIR filtering and multirate processing

## FIR convention

For coefficients \(h[0],\ldots,h[M-1]\), the causal FIR output is

\[
y[n] = \sum_{k=0}^{M-1} h[k]x[n-k],
\]

with samples before the start of a finite input treated as zero.

`fir::direct`, `fir::fft`, and `fir::overlap_save` return one causal output sample for every supplied input sample, so they omit the final \(M-1\) zero-input tail samples. `fir::full_direct` returns the complete finite linear convolution of lengths \(N\) and \(M\), with length \(N+M-1\).

`fir::Filter<T>` is the sample/block streaming direct form with a circular delay line. `fir::OverlapSaveFilter<T>` is the reusable FFT block form with an explicit transform size and precomputed kernel spectrum. Neither performs method selection.

General FIR execution supports binary32/binary64 real and complex samples. The symmetry-reduced methods are real linear-phase mechanisms.

## Linear-phase reductions

For a symmetric impulse response,

\[
h[k]=h[M-1-k],
\]

paired products become

\[
h[k]\bigl(x[n-k]+x[n-(M-1-k)]\bigr).
\]

`fir::symmetric_direct` evaluates this identity explicitly. For antisymmetric coefficients,

\[
h[k]=-h[M-1-k],
\]

`fir::antisymmetric_direct` evaluates

\[
h[k]\bigl(x[n-k]-x[n-(M-1-k)]\bigr).
\]

Both methods verify their structural precondition and fail instead of falling back to the ordinary direct sum.

## FFT FIR execution

`fir::fft` uses zero-padded radix-2 convolution. `fir::overlap_save` uses the explicit overlap-save implementation and requires the caller to choose the power-of-two transform size. Both preserve the causal FIR convention by retaining the first \(N\) convolution samples.

No crossover threshold, autotuning, or direct/FFT dispatcher is part of the FIR API.

## Windowed-sinc design

Frequencies in `fir_design` use cycles per sample, so Nyquist is \(0.5\). The ideal low-pass prototype is

\[
h[n]=2f_c\,\operatorname{sinc}\left(2f_c(n-\alpha)\right),
\qquad
\alpha=\frac{M-1}{2},
\]

with

\[
\operatorname{sinc}(x)=\frac{\sin(\pi x)}{\pi x}.
\]

The finite prototype is multiplied pointwise by the selected window. `lowpass_windowed_sinc` has a convenience overload using a symmetric Hann window. High-pass and band-stop designs use spectral inversion and therefore require an odd tap count; band-pass design subtracts two low-pass prototypes.

Design functions only produce coefficients. They do not construct or choose execution methods.

## Parks-McClellan / Remez exchange

`fir_design::remez_type1` implements odd-length Type-I linear-phase equiripple design. For \(M=2K+1\), the zero-phase response is

\[
A(f)=a_0+\sum_{k=1}^{K}a_k\cos(2\pi fk).
\]

At each exchange iteration, \(K+2\) extremal frequencies satisfy

\[
A(f_i)+\frac{(-1)^i\delta}{W(f_i)}=D(f_i),
\]

where \(D\) is the desired piecewise-constant response and \(W\) is the positive band weight. A dense grid evaluates the weighted error, alternating extrema are selected, and the system is solved again.

Transition regions are not approximation samples. Each specified band edge is treated as its own candidate extremum, so the exchange preserves the alternation problem across discontinuous pass/stop specifications.

The symmetric taps are reconstructed from the cosine coefficients as

\[
h[K]=a_0,\qquad h[K-k]=h[K+k]=\frac{a_k}{2}.
\]

`remez_lowpass` is the two-band convenience form.

## Elementary rate changes

`resampling::downsample(x,M,phase)` returns

\[
y[r]=x[\text{phase}+rM].
\]

`resampling::upsample_zero(x,L)` inserts \(L-1\) zeros between adjacent samples. A finite length-\(N\) input uses the minimal zero-stuffed representation

\[
(N-1)L+1,
\]

with no trailing zero group after the final original sample.

These primitives intentionally do not add anti-alias or anti-imaging filtering.

## FIR interpolation, decimation, and rational conversion

`interpolate_fir` performs zero insertion followed by the complete finite FIR convolution. `decimate_fir` performs the complete finite FIR convolution first and then downsamples. Keeping the complete finite response makes both operations algebraically consistent with the rational and polyphase forms.

The library does not multiply interpolation coefficients by \(L\) automatically; any passband-gain normalization belongs to the supplied coefficients.

For interpolation factor \(L\), decimation factor \(M\), and FIR \(h\), `rational_fir` performs exactly:

1. zero-stuff by \(L\);
2. linearly convolve with \(h\);
3. keep every \(M\)-th output sample.

## Polyphase rational conversion

Writing

\[
h_p[q]=h[p+qL],\qquad 0\le p<L,
\]

and considering output sample \(r\) at high-rate index \(n=rM\), define

\[
p=n\bmod L,\qquad c=\left\lfloor\frac{n}{L}\right\rfloor.
\]

Then the zero terms of explicit upsampling can be removed algebraically:

\[
y[r]=\sum_q h_p[q]x[c-q].
\]

`polyphase_rational` evaluates this expression directly, so it is equivalent to `rational_fir` without materializing or multiplying by inserted zeros. `polyphase_decompose` exposes the phase decomposition, while `polyphase_interpolate` and `polyphase_decimate` are the \(M=1\) and \(L=1\) special cases.

Rate factors and coefficients remain explicit. The API does not silently reduce rational factors, design a filter, or choose another implementation.

## Elementary interpolation algorithms

The repository also exposes simpler interpolation mechanisms separately:

- `zero_order_hold` repeats each sample;
- `linear` inserts points on straight-line segments;
- `windowed_sinc` evaluates finite Lanczos-windowed Whittaker-Shannon reconstruction with an explicit sinc radius.

The finite sinc form preserves original sample positions up to floating-point rounding but remains an approximation to infinite ideal sinc reconstruction.

## Numerical behavior

Execution and resampling use binary32/binary64 real or complex arithmetic. Different FIR organizations change operation order and therefore need not be bit-identical. Correctness tests compare direct, FFT, overlap-save, streaming, and polyphase forms with format-specific tolerances.

The Remez system uses partial pivoting and evaluates the exchange grid in the requested floating-point format. Very high-order or ill-conditioned specifications may therefore fail explicitly rather than being replaced by a different design method.
