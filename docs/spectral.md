# Time-frequency and spectral analysis

This repository keeps transform coefficients, spectral estimators, and time-frequency analysis as separate algorithms. Functions that depend on an FFT mechanism name the mechanism explicitly; the current arbitrary-length implementations use Bluestein.

## STFT convention

For a frame length \(L\), hop \(H\), window \(w[n]\), and frame index \(m\), the short-time Fourier transform is

\[
X_m[k] = \sum_{n=0}^{L-1} x[mH+n]w[n]e^{-2\pi i kn/L}.
\]

`stft::bluestein` evaluates every frame with the repository's arbitrary-length Bluestein FFT. The returned frames contain all \(L\) complex DFT bins. The public hop and window are never selected automatically.

With `pad_end=false`, only complete frames are emitted. With `pad_end=true`, frame starts continue while the start offset is inside the finite input and missing samples in the last frames are zero. No centering or implicit leading padding is performed.

Real and complex binary32/binary64 sample streams are supported. The analysis window is always real-valued in the corresponding scalar format.

## Inverse STFT and weighted overlap-add

`stft::inverse_bluestein` transforms every spectrum back to the time domain and performs weighted overlap-add with the supplied window. With the same analysis and synthesis window, the reconstruction is

\[
\hat{x}[n] =
\frac{\sum_m w[n-mH]\,\operatorname{IDFT}(X_m)[n-mH]}
     {\sum_m w^2[n-mH]}.
\]

The denominator is accumulated sample by sample rather than assuming a constant-overlap-add window. This makes the normalization explicit and works for any window/hop combination wherever the denominator is nonzero.

If a requested output sample has zero overlap weight, the result at that sample is zero. The library does not invent boundary padding or another synthesis window to make an otherwise uncovered sample reconstructible. Callers that require perfect finite-boundary reconstruction must choose padding/window/hop conventions that give every requested sample nonzero overlap weight.

`stft::power_spectrogram` is only the raw squared magnitude \(|X_m[k]|^2\) of already-computed STFT frames. Its one-sided option selects the unique nonnegative-frequency bins of real input; it does not apply PSD density scaling or one-sided energy doubling. The density-scaled spectrogram is a separate estimator in `spectral`.

## Magnitude, phase, and power spectra

For an existing complex spectrum \(X[k]\):

\[
M[k] = |X[k]|,\qquad
\phi[k] = \arg X[k],\qquad
P[k] = |X[k]|^2.
\]

`spectral::magnitude`, `spectral::phase`, and `spectral::power` apply exactly these definitions. They perform no FFT, windowing, density normalization, decibel conversion, or one-sided adjustment.

## Periodogram power spectral density

For a windowed length-\(L\) segment and sample rate \(F_s\), define

\[
U = \sum_{n=0}^{L-1} w^2[n].
\]

`spectral::periodogram_bluestein` uses

\[
S_{xx}[k] = \frac{|X[k]|^2}{F_s U}.
\]

This is power spectral density scaling, with frequency-bin spacing

\[
\Delta f = \frac{F_s}{L}.
\]

For a rectangular window, summing a two-sided periodogram and multiplying by \(\Delta f\) recovers the finite-signal mean-square value through Parseval's identity.

For real input, `Sides::one_sided` returns bins \(0\) through \(\lfloor L/2\rfloor\). Interior positive-frequency bins are doubled so that integrating the one-sided density preserves the same total power. DC is never doubled, and the Nyquist bin is not doubled when \(L\) is even. One-sided density is rejected for complex input because positive and negative frequencies are not redundant.

The window is explicit and must have positive energy. No default window or detrending operation is hidden in the estimator.

## Welch PSD

Welch's method averages periodograms of overlapping segments:

\[
\hat S_{xx}[k] = \frac{1}{R}\sum_{r=0}^{R-1}S^{(r)}_{xx}[k].
\]

`spectral::welch_bluestein` uses the caller-supplied window and hop. Only complete frames are included; trailing partial data is not silently padded. If the input is shorter than one window, the result is empty.

The estimator deliberately does not remove means, detrend, choose overlap, or select a window automatically. Those choices materially change the estimator and remain visible to the caller.

## Cross-spectral density

For two equally sized signals, the per-frame cross spectrum is

\[
S_{xy}^{(r)}[k] =
\frac{X_r[k]\overline{Y_r[k]}}{F_s U}.
\]

`spectral::cross_spectral_density_bluestein` averages this complex quantity over complete frames using the same explicit hop/window rules as Welch PSD. Setting \(x=y\) therefore yields the corresponding Welch power spectral density, up to floating-point rounding.

The complex result retains relative phase information. For real signals the same one-sided doubling rules used by the PSD estimator are applied to the nonnegative-frequency bins.

## Spectrogram density

`spectral::spectrogram_bluestein` applies the periodogram density normalization independently to every STFT frame, producing a time-frequency PSD matrix. Unlike Welch, it does not average frames.

`pad_end` is explicit. `false` keeps only complete frames; `true` uses the STFT zero-padded final-frame convention. One-sided real spectrograms use the same DC/Nyquist/interior-bin scaling as the one-sided periodogram.

## Goertzel algorithm

For angular frequency \(\omega=2\pi f\), Goertzel's second-order recurrence is

\[
s_n = x_n + 2\cos(\omega)s_{n-1} - s_{n-2},
\qquad s_{-1}=s_{-2}=0.
\]

Many Goertzel implementations need only \(|X|^2\), so they do not expose the absolute phase convention of the final state. This repository returns the complex Fourier coefficient itself. After processing \(N\) samples,

\[
X(f) = e^{-i\omega(N-1)}s_{N-1} - e^{-i\omega N}s_{N-2}.
\]

The terminal rotations are therefore part of `goertzel::frequency`; they make the returned complex coefficient agree with the repository's forward DFT convention rather than merely matching its magnitude.

`goertzel::bin(x,k)` uses \(f=k/N\) and is cross-checked against the corresponding DFT bin. `goertzel::frequency` accepts any explicit frequency in \([0,1)\) cycles per sample. Both support real and complex binary32/binary64 samples.

## Hilbert transform and analytic signal

For a real finite sequence, the analytic signal is formed in the DFT domain by suppressing negative frequencies and doubling the positive frequencies. Let \(X[k]\) be the length-\(N\) DFT. The multiplier is

- \(1\) at DC;
- \(2\) at strictly positive frequencies below Nyquist;
- \(1\) at the Nyquist bin when \(N\) is even;
- \(0\) at negative-frequency bins.

Then

\[
z[n] = \operatorname{IDFT}\{H[k]X[k]\}
     = x[n] + i\,\mathcal{H}\{x\}[n].
\]

`hilbert::analytic_bluestein` returns \(z[n]\). `hilbert::transform_bluestein` returns its imaginary component, and `hilbert::envelope_bluestein` returns \(|z[n]|\).

These are finite-block DFT-domain Hilbert transforms. Their boundary behavior is therefore periodic with the finite transform block; the library does not silently extend, reflect, or pad the input.

## Numerical behavior

All public algorithms in this milestone use binary32 or binary64. FFT-based methods inherit the floating-point behavior and arbitrary-length Bluestein contract documented for the Fourier-transform catalog.

STFT reconstruction, Welch averaging, and cross-spectral averaging change floating-point operation order relative to a single direct expression, so correctness is based on format-specific tolerances and mathematical identities rather than bit equality. The tests cover real and complex STFT reconstruction, PSD/CSD consistency, one-sided power conservation for a bin-centered real sinusoid, Goertzel phase agreement with the DFT, arbitrary-frequency Goertzel evaluation, and analytic-signal reconstruction for a bin-centered cosine.

The repository intentionally does not include timing, estimator campaign statistics, automatic window selection, automatic overlap selection, or spectral plotting. Those belong outside the algorithm library.
