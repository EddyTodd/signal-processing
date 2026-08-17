# Mathematics

This repository keeps the mathematical definition of each operation separate from its implementation. Method names in the public API identify the algorithm actually executed.

## Discrete Fourier transform

For a length-\(N\) sequence \(x_n\),

\[
X_k = \sum_{n=0}^{N-1} x_n e^{-2\pi i kn/N}.
\]

The inverse uses the positive exponential and the normalization \(1/N\). `fft::dft` evaluates this definition directly in \(O(N^2)\) and serves as the small-size correctness reference.

### Power-of-two mechanisms

`fft::radix2` is the canonical iterative radix-2 Cooley-Tukey decomposition with bit-reversal followed by butterfly stages. `fft::radix2_recursive` expresses the same factorization recursively by separating even and odd samples. `fft::stockham` is the radix-2 autosort formulation: each stage writes to a second buffer and therefore avoids a separate bit-reversal permutation.

`fft::radix4` decomposes the transform into four residue classes at each recursive level. `fft::split_radix` combines one \(N/2\) transform with two \(N/4\) transforms. `fft::modified_split_radix` implements the Johnson-Frigo scaled conjugate-pair mechanism: recursive scale factors turn selected complex twiddle products into tangent/cotangent forms with fewer real multiplications. It is retained as a faithful structural implementation rather than a claim that readable recursive C++ reaches the instruction count of generated machine-specific codelets.

### Composite and arbitrary lengths

For \(N=rm\), mixed-radix Cooley-Tukey recursively transforms \(m\)-point subsequences and combines them with an \(r\)-point DFT and cross-stage twiddle factors. `fft::mixed_radix` chooses the smallest available factor recursively.

When \(N=ab\) and \(\gcd(a,b)=1\), Good-Thomas / the Prime Factor Algorithm uses the Chinese Remainder Theorem to replace the cross-stage twiddles with input/output permutations. `fft::good_thomas` therefore performs independent \(a\)- and \(b\)-point transforms with no top-level twiddle multiplication.

For prime \(N\), Rader maps the nonzero indices to the multiplicative group modulo \(N\), reducing the DFT to a cyclic convolution of length \(N-1\). `fft::rader` performs that reduction explicitly.

Bluestein uses

\[
kn = \frac{k^2+n^2-(k-n)^2}{2}
\]

to rewrite the DFT as a convolution of quadratic chirps. `fft::bluestein` therefore supports any transform length. The chirp \(e^{\pm i\pi k^2/N}\) is periodic in \(k^2\) modulo \(2N\). The implementation computes that residue with overflow-safe modular integer multiplication before converting to `float` or `double`; it does not use runtime `long double`, whose width differs across common ABIs.

### Codelets and reusable plans

`SmallDftCodelet<T>` supports radices 2, 3, 4, 5, and 7. Radix-2 and radix-4 use explicit butterflies. Radix-3/5/7 precompute their root matrices so repeated execution performs no trigonometric setup.

The reusable plans are algorithm-specific:

- `Radix2Plan<T>` stores bit-reversal indices and twiddles;
- `RealRadix2Plan<T>` packs an even-length real signal into a half-size complex FFT and reconstructs the \(N/2+1\) unique bins implied by Hermitian symmetry;
- `MixedRadixPlan<T>` stores the factor tree, stage twiddles, and codelets;
- `GoodThomasPlan<T>` stores the CRT permutations and exposes zero top-level twiddles;
- `RaderPlan<T>` stores the prime-index permutations and transformed convolution kernel;
- `BluesteinPlan<T>` stores the chirp and transformed convolution kernel.

Plan construction may allocate and generate tables. Execution uses caller-provided scratch where needed and performs no hidden algorithm selection. This distinction is algorithmic: decomposition tables, permutations, chirps, and scratch requirements are part of a reusable FFT implementation rather than benchmark infrastructure.

### Explicit machine kernels

`KernelRadix2Plan` is an explicit binary64 extension with scalar, AVX2/FMA, and AVX-512/FMA kernels. There is no automatic fallback. A requested vector kernel is constructed only when that implementation is compiled for the target and the processor reports the required ISA; otherwise construction fails. This keeps ISA identity inspectable for external benchmarking.

### FFT numerical conventions

Only IEEE-oriented binary32 (`float`) and binary64 (`double`) are public FFT scalar types. Forward transforms are unnormalized. Every inverse transform applies exactly one \(1/N\) normalization. Cross-method tests compare the structural algorithms with the direct DFT and separately check forward/inverse reconstruction for both formats.

The implementations accumulate ordinary floating-point error through additions, twiddle multiplication, cancellation, and stage composition, so the tests use format-specific error budgets rather than requiring bit identity. Architecture-specific FMA kernels can legitimately round differently from scalar complex multiplication while still satisfying the same transform contract.

## Convolution

For sequences \(x\) and \(h\), linear convolution is

\[
y[n] = (x*h)[n] = \sum_k x[k]h[n-k].
\]

There is no conjugation in convolution, so the same definition applies to real and complex samples. For finite inputs of lengths \(N\) and \(M\), the full linear result has length \(N+M-1\).

`convolution::direct` evaluates this definition in \(O(NM)\) and is the correctness reference. `convolution::fft` chooses a power-of-two transform length \(K\ge N+M-1\), zero-pads both operands, computes

\[
Y = \operatorname{DFT}_K(x)\,\operatorname{DFT}_K(h),
\]

pointwise, and applies the inverse transform. The padding is what converts the FFT's circular convolution into the desired linear convolution.

For equal-length sequences, circular convolution is

\[
y[n] = \sum_{k=0}^{N-1}x[k]h[(n-k)\bmod N].
\]

`convolution::circular` evaluates the periodic sum directly. `convolution::circular_bluestein` transforms exactly \(N\) samples with Bluestein, so arbitrary circular lengths are represented without padding to a different circular period.

### Overlap-add

Let the input be split into blocks of \(L\) new samples,

\[
x[n] = \sum_b x_b[n-bL].
\]

Linearity gives

\[
x*h = \sum_b (x_b*h)[n-bL].
\]

`convolution::overlap_add` FFT-convolves each block with the fixed kernel using a transform large enough for that block's complete linear result, then adds the overlapping tails at their shifted positions. The public `block_size` parameter is part of the method; the library does not auto-tune it.

### Overlap-save

For a kernel of length \(M\) and an FFT length \(K\ge M\), overlap-save consumes

\[
L = K-M+1
\]

new samples per block. Each transform input consists of the previous \(M-1\) time samples followed by those \(L\) new samples. The length-\(K\) circular convolution contains aliasing in its first \(M-1\) outputs, so those outputs are discarded and the remaining \(L\) are valid linear-convolution samples.

`OverlapSavePlan<T>` stores the fixed kernel spectrum and the \(M-1\)-sample input overlap. `convolution::overlap_save` pads the stream with initial and final zeros as needed to return the complete \(N+M-1\) linear result. The transform length is explicit and must be a power of two because this implementation deliberately uses `Radix2Plan`.

### Streaming direct convolution

`StreamingDirect<T>` keeps the last \(M\) input samples in a circular history buffer and evaluates the defining convolution sum once per new sample. `process()` therefore emits the causal output for the newly supplied sample or block; `flush()` feeds the \(M-1\) trailing zeros required to recover the tail of a finite linear convolution.

### Uniform partitioned convolution

For long fixed kernels, write the kernel as equal time-domain partitions of length \(P\):

\[
h[n] = \sum_{p=0}^{Q-1} h_p[n-pP].
\]

Input is processed in matching \(P\)-sample blocks. If \(X_b\) is the FFT of input block \(b\) and \(H_p\) is the precomputed FFT of kernel partition \(p\), a frequency-domain delay line forms

\[
Y_b[k] = \sum_{p=0}^{Q-1} X_{b-p}[k]H_p[k].
\]

One inverse FFT produces the current block contribution, and its second half is overlap-added into the next block. `UniformPartitionedConvolver<T>` stores the kernel spectra, input-spectrum history, and time-domain overlap; `convolution::partitioned` is the finite-signal wrapper. `partition_size` remains explicit because latency, memory, and arithmetic work are properties of this algorithmic choice rather than library policy.

### Convolution numerical behavior

All convolution forms support binary32/binary64 real samples and `std::complex<float>` / `std::complex<double>`. FFT-based real convolution performs the arithmetic in the corresponding complex format and returns the real component after inversion; the discarded imaginary component is roundoff, not a second signal channel. Different blocking or partitioning choices change the order of floating-point operations, so cross-method tests use format-specific tolerances instead of requiring bit identity.

## Correlation

The library uses the lag convention

\[
r_{xy}[\ell] = \sum_n x[n]\overline{y[n-\ell]}.
\]

The complex conjugate is essential for complex signals and disappears for real signals. The returned vector stores lag \(\ell\) at index

\[
\ell + (|y|-1),
\]

so zero lag is at index \(|y|-1\).

`correlation::cross_direct` evaluates the sum directly. Reversing and conjugating the second sequence gives

\[
r_{xy} = x * \operatorname{reverse}(\overline{y}),
\]

which is the identity used by `correlation::cross_fft`. `auto_direct` and `auto_fft` set \(y=x\); the initial `cross` and `auto_correlation` names remain compatibility aliases for the direct definitions.

For every lag, normalized cross-correlation divides by the energies of the samples that actually overlap at that lag:

\[
\rho_{xy}[\ell] =
\frac{r_{xy}[\ell]}
{\sqrt{\left(\sum_{n\in I_\ell}|x[n]|^2\right)
       \left(\sum_{n\in I_\ell}|y[n-\ell]|^2\right)}}.
\]

`normalized_direct` uses the direct numerator; `normalized_fft` uses the FFT numerator. Both compute the overlap energies with prefix sums. If either overlapping segment has zero energy, the normalized value is defined as zero. This is an energy-normalized correlation coefficient, not the separate biased/unbiased lag-count normalization used by some statistical autocorrelation conventions.

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

The library is an algorithm library rather than a benchmark framework. No method silently substitutes another implementation. Stateful objects own only state intrinsic to their algorithm; benchmarking, timing, statistics, campaign orchestration, and reporting remain outside this repository.
