# IIR filtering and filter design

This document fixes the transfer-function, realization, prototype, and frequency conventions used by the IIR API. Execution structures and design algorithms remain separate public mechanisms.

## Difference equation

For numerator coefficients \(b_0,\ldots,b_M\) and denominator coefficients \(a_0,\ldots,a_N\), the library normalizes \(a_0=1\) and implements

\[
y[n]=\sum_{k=0}^{M} b_k x[n-k]-\sum_{k=1}^{N} a_k y[n-k].
\]

Equivalently,

\[
H(z)=\frac{\sum_{k=0}^{M} b_k z^{-k}}
           {1+\sum_{k=1}^{N} a_k z^{-k}}.
\]

The general `DirectFormI`, `DirectFormII`, `TransposedDirectFormI`, and `TransposedDirectFormII` classes implement the same transfer function with different state graphs. No class silently changes realization.

- Direct Form I keeps separate input and output delay lines.
- Direct Form II uses one canonical delay line. With internal state \(w\),

  \[
  w[n]=x[n]-\sum_{k=1}^{N}a_k w[n-k],\qquad
  y[n]=\sum_{k=0}^{M}b_k w[n-k].
  \]

- The transposed classes are the signal-flow transposes of their corresponding direct structures.

See [Infinite impulse response](https://en.wikipedia.org/wiki/Infinite_impulse_response).

## First- and second-order sections

A biquad has

\[
H(z)=\frac{b_0+b_1z^{-1}+b_2z^{-2}}
           {a_0+a_1z^{-1}+a_2z^{-2}}.
\]

`BiquadDirectFormI`, `BiquadDirectFormII`, and `BiquadTransposedDirectFormII` expose the three realizations explicitly. `FirstOrderTransposedDirectFormII` represents a first-order section without padding it conceptually into a second-order filter.

For a cascade of \(Q\) sections,

\[
H(z)=g\prod_{q=0}^{Q-1} H_q(z).
\]

`SosCascade` executes that expression as transposed Direct Form II biquads. The overall gain \(g\) remains separate from the section coefficients. See [Digital biquad filter](https://en.wikipedia.org/wiki/Digital_biquad_filter).

For high-order filters, an SOS representation is generally the numerically preferable execution form because it avoids forming one high-degree denominator polynomial. The full transfer-function coefficients remain available for study and cross-checking; the library does not hide either representation.

## Pole-zero-gain representation

Filter design uses explicit pole-zero-gain (ZPK) form. An analog prototype is

\[
H(s)=k\frac{\prod_i(s-z_i)}{\prod_j(s-p_j)},
\]

and a digital filter is represented by the corresponding roots in the \(z\)-plane. Analog stability requires

\[
\operatorname{Re}(p_j)<0,
\]

while digital stability requires

\[
|p_j|<1.
\]

`analog_stable` and `stable` check those conditions directly.

## Analog low-pass prototypes

All prototype functions use a normalized low-pass analog form. The public family name is always explicit.

### Butterworth

`butterworth_prototype` places the poles uniformly on the left half of the unit circle, giving the maximally flat magnitude response. The normalized critical frequency is the \(-3\) dB point. See [Butterworth filter](https://en.wikipedia.org/wiki/Butterworth_filter).

### Chebyshev Type I

For passband ripple \(r_p\) dB,

\[
\epsilon=\sqrt{10^{r_p/10}-1}.
\]

`chebyshev1_prototype` maps the pole pattern through the hyperbolic ellipse determined by \(\epsilon\). See [Chebyshev filter](https://en.wikipedia.org/wiki/Chebyshev_filter).

### Chebyshev Type II

For stopband attenuation \(r_s\) dB, `chebyshev2_prototype` uses the inverse-Chebyshev pole construction and the finite imaginary-axis transmission zeros required by the Type-II response.

### Elliptic / Cauer

`elliptic_prototype` implements the Cauer construction from the passband and stopband ripple specifications. It solves the elliptic degree equation through the nome expansion and evaluates the required Jacobi quantities directly. The implementation uses Carlson's symmetric integral \(R_F\) internally for complete/incomplete elliptic integrals, then a monotone inverse for the real Jacobi amplitude so the inverse cannot jump to a different branch near modulus one.

The result has ripple in both passband and stopband and finite transmission zeros. See [Elliptic filter](https://en.wikipedia.org/wiki/Elliptic_filter).

### Bessel / Thomson

`bessel_prototype` obtains poles from the reverse Bessel polynomial using a scaled Aberth-Ehrlich simultaneous root iteration. Three normalizations are explicit:

- `phase`: phase-matched normalization, with Butterworth-like asymptotes;
- `delay`: the natural Bessel polynomial, with unit low-frequency group delay;
- `magnitude_3db`: scales the prototype so \(|H(j)|=1/\sqrt{2}\).

See [Bessel filter](https://en.wikipedia.org/wiki/Bessel_filter).

## Analog frequency transformations

A normalized low-pass prototype can be transformed explicitly before digitization. Let \(\Omega_0>0\) be a center/cutoff frequency and \(B>0\) a bandwidth.

Low-pass scaling uses

\[
s\mapsto \frac{s}{\Omega_0}.
\]

High-pass transformation uses

\[
s\mapsto \frac{\Omega_0}{s}.
\]

Band-pass transformation uses

\[
s\mapsto \frac{s^2+\Omega_0^2}{Bs},
\]

and band-stop uses

\[
s\mapsto \frac{Bs}{s^2+\Omega_0^2}.
\]

The ZPK functions preserve zeros at infinity through relative-degree bookkeeping instead of manufacturing finite roots prematurely.

## Bilinear transform and digital frequencies

Digital design uses the [bilinear transform](https://en.wikipedia.org/wiki/Bilinear_transform)

\[
s=2f_s\frac{z-1}{z+1}.
\]

The public digital design functions use normalized frequency \(f\in(0,1)\), where \(f=1\) is Nyquist. With the internal normalized sample rate \(f_s=1\), each requested digital critical frequency is prewarped as

\[
\Omega=2\tan\left(\frac{\pi f}{2}\right).
\]

The analog prototype is then transformed to low-pass, high-pass, band-pass, or band-stop form and mapped to the digital plane. Zeros at analog infinity become digital zeros at \(z=-1\).

`butterworth`, `chebyshev1`, `chebyshev2`, `elliptic`, and `bessel` expose this pipeline while retaining the family name in the call. `Response` selects only the frequency transformation; it is not an algorithm dispatcher.

## Coefficient and SOS conversion

`transfer_function` expands a balanced digital ZPK representation into real numerator/denominator coefficients. `second_order_sections` pairs conjugate real-coefficient roots into first/second-order factors and leaves the overall ZPK gain explicit in `SosDesign`.

These conversions deliberately reject a result whose imaginary residue is too large to be explained by binary32/binary64 roundoff. A malformed or non-conjugate root set therefore does not silently become a real filter.

## Numerical conventions

- IIR scalar APIs support binary32 and binary64 only.
- `a[0]` is normalized once at construction; a zero leading denominator coefficient is invalid.
- Analog prototype functions return finite zeros only; relative degree represents zeros at infinity.
- Digital critical frequencies are normalized to Nyquist, not sample rate.
- No filter family, realization, normalization, or SOS/transfer representation is chosen automatically.
- Bessel polynomial roots are solved after variable scaling to reduce high-order conditioning problems.
- Elliptic inverse-Jacobi work uses monotone inversion to remain on the intended real branch.
- Tests compare mathematically equivalent realizations and check family-defining response points rather than requiring bit-identical state histories.
