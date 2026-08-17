# FIR filtering and multirate processing

## FIR convention

For coefficients \(h[0],\ldots,h[M-1]\), the causal FIR output is

\[
y[n] = \sum_{k=0}^{M-1} h[k]x[n-k],
\]

with samples before the start of a finite input treated as zero.

`fir::direct`, `fir::fft`, and `fir::overlap_save` return one causal output sample for every supplied input sample. They therefore omit the final \(M-1\) zero-input tail samples. `fir::full_direct` returns the complete finite linear convolution of lengths \(N\) and \(M\), which has length \(N+M-1\).

`fir::Filter<T>` is the sample/block streaming direct form with a circular delay line. `fir::OverlapSaveFilter<T>` is the reusable block form whose transform size is explicit and whose kernel spectrum is precomputed. Neither class selects a method dynamically.

Real and complex binary32/binary64 samples are supported by the general FIR execution forms. The symmetric and antisymmetric reductions are real linear-phase mechanisms and therefore use real binary32/binary64 samples and coefficients.

## Linear-phase symmetry reductions

For a symmetric length-\(M\) impulse response,

\[
h[k]=h[M-1-k],
\]

paired products can be written

\[
h[k]\bigl(x[n-k]+x[n-(M-1-k)]\bigr).
\]

`fir::symmetric_direct` evaluates that pairing explicitly. For an antisymmetric response,

\[
h[k]=-h[M-1-k],
\]

`fir::antisymmetric_direct` instead evaluates

\[
h[k]\bigl(x[n-k]-x[n-(M-1-k)]\bigr).
\]

These methods verify their structural precondition and fail rather than silently falling back to the ordinary direct sum.

## FFT FIR execution

`fir::fft` uses the repository's explicit radix-2 zero-padded convolution mechanism. `fir::overlap_save` uses the explicit overlap-save implementation and requires the caller to choose the power-of-two transform size. These functions preserve the FIR causal-output convention by retaining the first \(N\) convolution samples.

No crossover threshold, autotuning, or hidden direct/FFT dispatch is part of the FIR API.

## Windowed-sinc FIR design

Frequencies in `fir_design` are normalized in cycles per sample, so Nyquist is \(0.5\). The ideal low-pass impulse response with cutoff \(f_c\) is sampled as

\[
h[n] = 2f_c\,\operatorname{sinc}\left(2f_c(n-\alpha)\right),
\qquad
\alpha=\frac{M-1}{2},
\]

where

\[
\operatorname{sinc}(x)=\frac{\sin(\pi x)}{\pi x}.
\]

The finite sequence is multiplied pointwise by the caller-supplied window. `lowpass_windowed_sinc` also has a convenience overload using the repository's symmetric Hann window.

High-pass and band-stop designs use spectral inversion and therefore require an odd tap count so the unit impulse lies on an integer center sample. Band-pass design is the difference of two low-pass prototypes.

Filter design is deliberately separate from filter execution: design functions produce coefficient vectors and do not construct or choose a runtime filtering method.

## Parks-McClellan / Remez exchange

`fir_design::remez_type1` implements the odd-length Type-I linear-phase exchange algorithm. For \(M=2K+1\), the zero-phase response is represented as

\[
A(f)=a_0+\sum_{k=1}^{K}a_k\cos(2\pi fk).
\]

At each exchange iteration, \(K+2\) extremal frequencies satisfy

\[
A(f_i)+\frac{(-1)^i\delta}{W(f_i)}=D(f_i),
\]

where \(D\) is the desired piecewise-constant response and \(W\) is the positive band weight. The dense-grid weighted error is then searched for alternating extrema and the system is solved again.

Transition regions between specified bands are not part of the approximation grid. Each pass/stop-band edge is considered independently when extrema are selected, which preserves the alternation problem across discontinuous desired responses.

The returned Type-I taps are reconstructed from the cosine coefficients as

\[
h[K]=a_0,\qquad h[K-k]=h[K+k]=\frac{a_k}{2}.
\]

`remez_lowpass` is a two-band convenience wrapper. The implementation is a readable exchange algorithm, not an external-library wrapper or an autotuned filter designer.

## Downsampling and zero insertion

`resampling::downsample(x, M, phase)` returns

\[
y[r]=x[\text{phase}+rM].
\]

`resampling::upsample_zero(x, L)` inserts \(L-1\) zeros between adjacent input samples. For a finite length-\(N\) input it uses the minimal zero-stuffed representation of length

\[
(N-1)L+1,
\]

so no artificial trailing zero group is appended after the final original sample.

These primitives do not include an anti-alias or anti-imaging filter. The filtered forms make that FIR operation explicit.

## FIR decimation and interpolation

`decimate_fir` first applies the causal FIR and then downsamples. `interpolate_fir` performs zero insertion followed by the complete FIR convolution. The library does not multiply interpolation coefficients by \(L\) automatically; any passband-gain normalization belongs to the supplied filter coefficients.

This choice keeps the operation algebraically transparent and makes direct and polyphase forms exactly comparable.

## Rational sample-rate conversion

For interpolation factor \(L\), decimation factor \(M\), and FIR \(h\), `rational_fir` performs the literal sequence

1. zero-stuff by \(L\);
2. linearly convolve with \(h\);
3. keep every \(M\)-th output sample.

`polyphase_rational` computes the same sequence without materializing or multiplying by inserted zeros. Writing

\[
h_p[q]=h[p+qL],\qquad 0\le p<L,
\]

and considering output sample \(r\) at high-rate index \(n=rM\), let

\[
p=n\bmod L,\qquad c=\left\lfloor\frac{n}{L}\right\rfloor.
\]

Then

\[
y[r]=\sum_q h_p[q]x[c-q],
\]

with out-of-range input indices omitted. `polyphase_decompose` exposes the phase decomposition directly; `polyphase_interpolate` and `polyphase_decimate` are the corresponding special cases.

The public factors and coefficient set are always explicit. There is no rational-factor reduction, filter synthesis, or algorithm selection hidden inside the conversion call.

## Elementary interpolation

The repository also keeps simple interpolation mechanisms as distinct algorithms:

- `zero_order_hold` repeats every sample;
- `linear` inserts points on straight line segments;
- `windowed_sinc` evaluates a finite Lanczos-windowed Whittaker-Shannon reconstruction with an explicit sinc half-width.

The windowed-sinc implementation preserves original sample values at integer sample positions up to floating-point rounding. Its finite support makes it an approximation to ideal infinite sinc reconstruction, and the half-width remains an explicit algorithm parameter.

## Numerical behavior

All execution and resampling paths use binary32/binary64 real or complex arithmetic. Different FIR organizations change operation order and therefore need not be bit-identical. Correctness tests compare direct, FFT, overlap-save, streaming, and polyphase forms with format-specific tolerances.

The Remez linear system uses partial pivoting and the exchange grid is evaluated in the requested floating-point format. Large/high-order or extremely ill-conditioned specifications can therefore fail explicitly rather than being silently replaced by a different design method.
