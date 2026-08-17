# Mathematics

This repository keeps the mathematical definition of each operation separate from its implementation. Method names in the public API identify the algorithm actually executed.

## Discrete Fourier transform

For a length-\(N\) sequence \(x_n\),

\[
X_k = \sum_{n=0}^{N-1} x_n e^{-2\pi i kn/N}.
\]

The inverse uses the positive exponential and the normalization \(1/N\). `fft::dft` evaluates this definition directly in \(O(N^2)\). `fft::radix2` is the iterative radix-2 Cooley-Tukey decomposition and therefore requires a power-of-two length.

## Convolution

Linear convolution is

\[
y_n = \sum_k x_k h_{n-k}.
\]

`convolution::direct` evaluates the sum directly. `convolution::fft` zero-pads both operands, transforms them, multiplies pointwise, and applies the inverse transform according to the convolution theorem.

Circular convolution of equal-length sequences is

\[
y_n = \sum_{k=0}^{N-1} x_k h_{(n-k)\bmod N}.
\]

## Correlation

Cross-correlation is implemented using the lag convention

\[
r_{xy}[\ell] = \sum_n x[n]y[n-\ell].
\]

Autocorrelation is the special case \(x=y\).

## Discrete cosine transforms

The initial catalog implements the unnormalized DCT-I, II, III, and IV definitions directly. In particular,

\[
X_k^{\mathrm{II}} = 2\sum_{n=0}^{N-1}x_n
\cos\left[\frac{\pi}{N}\left(n+\frac12\right)k\right],
\]

and

\[
X_k^{\mathrm{III}} = x_0 + 2\sum_{n=1}^{N-1}x_n
\cos\left[\frac{\pi}{N}n\left(k+\frac12\right)\right].
\]

With these conventions, applying DCT-III to DCT-II returns \(2N x\). Fast FFT-based forms belong in a later algorithm-specific extension; the direct forms remain as definitions and correctness references.

## FIR filters

For coefficients \(b_k\), an FIR filter computes

\[
y[n] = \sum_{k=0}^{M-1} b_k x[n-k].
\]

`fir::direct` processes a finite block with zero initial history. `fir::Filter` implements the same recurrence as a stateful circular delay line for streaming use.

## IIR filters

For feed-forward coefficients \(b_k\) and feedback coefficients \(a_k\), normalized so \(a_0=1\),

\[
y[n] = \sum_{k=0}^{M} b_k x[n-k] - \sum_{k=1}^{N} a_k y[n-k].
\]

`iir::DirectFormI` represents this recurrence directly. `iir::BiquadTransposedDirectFormII` implements a second-order section in transposed Direct Form II, a common building block for stable higher-order filter cascades.

## Windows

For \(0\le n<N\), the Hann window is

\[
w[n] = \frac12 - \frac12\cos\left(\frac{2\pi n}{N-1}\right).
\]

The initial catalog also includes Hamming and Blackman windows. The definitions return a single value of one when \(N=1\), avoiding the otherwise undefined division by \(N-1\).

## Numerical conventions

The initial library supports `float` and `double` through constrained templates. No algorithm silently substitutes another method. FFT inverse normalization is applied exactly once by the inverse transform. Filtering objects own only their algorithmic state; benchmarking, timing, statistics, and campaign infrastructure are deliberately outside this repository.
