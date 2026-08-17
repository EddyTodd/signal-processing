# Discrete transforms and windows

This document records the exact transform conventions used by the public APIs. The direct implementations are the mathematical references; optimized variants retain the same normalization and output ordering.

## DCT-I through DCT-IV

For a real sequence \(x_0,\ldots,x_{N-1}\), the library uses the unnormalized cosine-transform conventions

\[
C_k^{\mathrm I}=x_0+(-1)^k x_{N-1}+2\sum_{n=1}^{N-2}x_n\cos\frac{\pi nk}{N-1},
\]

\[
C_k^{\mathrm {II}}=2\sum_{n=0}^{N-1}x_n
\cos\left[\frac{\pi}{N}\left(n+\frac12\right)k\right],
\]

\[
C_k^{\mathrm {III}}=x_0+2\sum_{n=1}^{N-1}x_n
\cos\left[\frac{\pi}{N}n\left(k+\frac12\right)\right],
\]

and

\[
C_k^{\mathrm {IV}}=2\sum_{n=0}^{N-1}x_n
\cos\left[\frac{\pi}{N}\left(n+\frac12\right)\left(k+\frac12\right)\right].
\]

`dct::type1`, `type2`, `type3`, and `type4` evaluate these definitions directly. DCT-I is self-inverse up to \(2(N-1)\) for \(N>1\); DCT-II and DCT-III are inverse pairs up to \(2N\); DCT-IV is self-inverse up to \(2N\).

The `_bluestein` forms reduce each definition to a Fourier transform while naming Bluestein explicitly. DCT-I uses an even extension of length \(2(N-1)\). DCT-II uses the even half-sample extension

\[
y_n=x_n,\qquad y_{2N-1-n}=x_n,
\]

followed by a phase rotation of the Fourier coefficients. DCT-III reconstructs the Hermitian spectrum corresponding to that extension and applies an inverse Bluestein transform. DCT-IV is obtained from the odd-index DCT-II coefficients of a zero-extended length-\(2N\) sequence. These constructions support arbitrary \(N\); there is no hidden power-of-two substitution.

## DST-I through DST-IV

The sine-transform conventions are

\[
S_k^{\mathrm I}=2\sum_{n=0}^{N-1}x_n
\sin\frac{\pi(n+1)(k+1)}{N+1},
\]

\[
S_k^{\mathrm {II}}=2\sum_{n=0}^{N-1}x_n
\sin\left[\frac{\pi}{N}\left(n+\frac12\right)(k+1)\right],
\]

\[
S_k^{\mathrm {III}}=(-1)^k x_{N-1}+
2\sum_{n=0}^{N-2}x_n
\sin\left[\frac{\pi}{N}(n+1)\left(k+\frac12\right)\right],
\]

and

\[
S_k^{\mathrm {IV}}=2\sum_{n=0}^{N-1}x_n
\sin\left[\frac{\pi}{N}\left(n+\frac12\right)\left(k+\frac12\right)\right].
\]

DST-I is self-inverse up to \(2(N+1)\); DST-II and DST-III are inverse pairs up to \(2N\); DST-IV is self-inverse up to \(2N\).

The `_bluestein` variants use odd extensions. DST-I places the input between two zero endpoints in a length-\(2(N+1)\) odd extension. DST-II uses a length-\(2N\) half-sample odd extension and a phase rotation. DST-III reconstructs the corresponding Hermitian spectrum, and DST-IV is obtained from the even-index DST-II coefficients of a zero-extended sequence.

## Discrete Hartley transform

The discrete Hartley transform is

\[
H_k=\sum_{n=0}^{N-1}x_n\operatorname{cas}\left(\frac{2\pi nk}{N}\right),
\qquad \operatorname{cas}(\theta)=\cos\theta+\sin\theta.
\]

It is real-to-real and self-inverse up to the factor \(N\). `hartley::direct` evaluates the definition in \(O(N^2)\). With the repository's forward DFT convention

\[
F_k=\sum_n x_n\left(\cos\theta-i\sin\theta\right),
\]

the identity

\[
H_k=\Re F_k-\Im F_k
\]

is used by `hartley::bluestein` for arbitrary lengths.

## Walsh-Hadamard transform

The library uses the natural Sylvester/Hadamard ordering rather than a hidden sequency permutation. For power-of-two \(N\),

\[
W_k=\sum_{n=0}^{N-1}(-1)^{\operatorname{popcount}(n\mathbin{\&}k)}x_n.
\]

`walsh_hadamard::direct` evaluates this matrix definition in \(O(N^2)\). `walsh_hadamard::fast` applies the standard butterfly

\[
(a,b)\mapsto(a+b,a-b)
\]

at successively doubled widths, giving \(O(N\log N)\) work. Applying either transform twice returns \(Nx\).

## MDCT and IMDCT

For an input block of length \(2N\), the modified discrete cosine transform is

\[
X_k=\sum_{n=0}^{2N-1}x_n
\cos\left[\frac{\pi}{N}\left(n+\frac12+\frac N2\right)
\left(k+\frac12\right)\right],\qquad 0\le k<N.
\]

The inverse definition used by `mdct::inverse_direct` is

\[
\tilde x_n=\frac{2}{N}\sum_{k=0}^{N-1}X_k
\cos\left[\frac{\pi}{N}\left(n+\frac12+\frac N2\right)
\left(k+\frac12\right)\right].
\]

A single MDCT block is critically sampled and therefore not independently invertible: the IMDCT contains time-domain aliasing. Perfect reconstruction is obtained in the usual lapped-transform setting with overlapping blocks and a compatible analysis/synthesis window. The library does not pretend that `inverse_direct(direct(x))` must equal one isolated input block.

## Window sampling convention

Most windows expose `windows::Sampling::symmetric` and `windows::Sampling::periodic` explicitly. Symmetric windows use denominator \(N-1\) and are appropriate when the endpoint symmetry itself matters, such as FIR design. Periodic windows use denominator \(N\), equivalent to taking the first \(N\) samples of a symmetric \(N+1\)-point cosine-series window, and are appropriate for an \(N\)-point periodic Fourier analysis.

The catalog contains rectangular, Bartlett, Hann, Hamming, Blackman, four-term Blackman-Harris, continuous-first-derivative Nuttall, flat-top, Kaiser, Gaussian, Tukey, Lanczos, and Welch windows.

The Kaiser window is

\[
w_n=\frac{I_0\!\left(\beta\sqrt{1-r_n^2}\right)}{I_0(\beta)},
\qquad r_n\in[-1,1],
\]

where \(I_0\) is evaluated by its convergent power series in the header. `beta=0` therefore reduces exactly to the rectangular window.

The Gaussian window uses

\[
w_n=\exp\left[-\frac12\left(\frac{r_n}{\sigma}\right)^2\right],
\]

where `sigma` is expressed relative to the normalized half-width \(r_n\in[-1,1]\). Tukey requires \(0\le\alpha\le1\); \(\alpha=0\) is rectangular and \(\alpha=1\) is Hann.

All window functions return an empty vector for \(N=0\) and the single value one for \(N=1\). No function silently changes the requested sampling convention or parameter values.
