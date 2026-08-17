# Mathematics

This repository keeps the mathematical definition of each operation separate from its implementation. Public method names identify the algorithm actually executed. Cross-library behavioral rules are frozen in [`contracts.md`](contracts.md).

## Discrete Fourier transform

For a length-\(N\) sequence \(x_n\),

\[
X_k = \sum_{n=0}^{N-1} x_n e^{-2\pi i kn/N}.
\]

The inverse uses the positive exponential and normalization \(1/N\). `fft::dft` evaluates the definition directly in \(O(N^2)\) and serves as the small-size correctness reference.

### Power-of-two mechanisms

`fft::radix2` is iterative radix-2 Cooley-Tukey with bit reversal followed by butterfly stages. `fft::radix2_recursive` expresses the same factorization recursively by separating even and odd samples. `fft::stockham` is the radix-2 autosort formulation and avoids a separate bit-reversal permutation by alternating source and destination layouts at each stage.

`fft::radix4` decomposes into four residue classes. `fft::split_radix` combines one \(N/2\) transform with two \(N/4\) transforms. `fft::modified_split_radix` implements the Johnson-Frigo scaled conjugate-pair mechanism: recursive scale factors convert selected complex twiddle products into tangent/cotangent forms with fewer real multiplications. It is retained as a faithful structural implementation rather than a claim that readable recursive C++ reaches the instruction count of generated machine-specific codelets.

### Composite and arbitrary lengths

For \(N=rm\), mixed-radix Cooley-Tukey recursively transforms \(m\)-point subsequences and combines them with an \(r\)-point DFT and cross-stage twiddle factors. `fft::mixed_radix` chooses the smallest available factor recursively.

When \(N=ab\) and \(\gcd(a,b)=1\), Good-Thomas / the Prime Factor Algorithm uses the Chinese Remainder Theorem to replace cross-stage twiddles with permutations. If \(b^{-1}\) denotes the inverse of \(b\) modulo \(a\), and \(a^{-1}\) the inverse of \(a\) modulo \(b\), one CRT input index can be written

\[
n \equiv n_1 b(b^{-1}\!\!\pmod a)+n_2 a(a^{-1}\!\!\pmod b) \pmod N.
\]

`fft::good_thomas(input, direction, a, b)` evaluates the resulting independent \(a\)- and \(b\)-point transforms with no top-level twiddle multiplication. The modular products and sums are performed without overflowing before reduction.

For prime \(N\), Rader maps nonzero indices to the multiplicative group modulo \(N\), reducing the DFT to cyclic convolution of length \(N-1\). `fft::rader` performs that reduction explicitly.

Bluestein uses

\[
kn = \frac{k^2+n^2-(k-n)^2}{2}
\]

to rewrite the DFT as convolution of quadratic chirps. `fft::bluestein` therefore supports every positive transform length. The chirp \(e^{\pm i\pi k^2/N}\) is periodic in \(k^2\) modulo \(2N\). The implementation computes that residue with overflow-safe modular integer multiplication before conversion to `float` or `double`; it does not use runtime `long double`.

### Codelets and reusable plans

`SmallDftCodelet<T>` supports radices 2, 3, 4, 5, and 7. Radix-2 and radix-4 use explicit butterflies. Radix-3/5/7 precompute root matrices so repeated execution performs no trigonometric setup.

Reusable plans are algorithm-specific:

- `Radix2Plan<T>` stores bit-reversal indices and twiddles;
- `RealRadix2Plan<T>` packs an even-length real signal into a half-size complex FFT and reconstructs the \(N/2+1\) unique bins implied by Hermitian symmetry;
- `MixedRadixPlan<T>` stores the factor tree, stage twiddles, and codelets;
- `GoodThomasPlan<T>` stores the CRT permutations and zero top-level twiddles;
- `RaderPlan<T>` stores prime-index permutations and the transformed convolution kernel;
- `BluesteinPlan<T>` stores the chirp and transformed convolution kernel.

Plan construction may allocate and generate tables. Execution uses caller-provided scratch where needed and performs no hidden algorithm selection. Required logical table and scratch sizes are checked before wrapped `std::size_t` arithmetic can become an allocation or index.

### Explicit machine kernels

`KernelRadix2Plan` is an explicit binary64 extension with scalar, AVX2/FMA, and AVX-512/FMA kernels. There is no automatic fallback. A requested vector kernel is constructed only when that implementation is compiled for the target and the processor reports the required ISA; otherwise construction fails.

### FFT numerical conventions

Only binary32 (`float`) and binary64 (`double`) are public FFT scalar types. Forward transforms are unnormalized. Every inverse FFT-family transform applies exactly one \(1/N\) normalization.

Different faithful decompositions accumulate floating-point error in different orders through additions, twiddle multiplication, cancellation, and stage composition. Cross-method tests therefore compare against the direct DFT with format-specific error budgets rather than requiring bit identity.

## Convolution

For sequences \(x\) and \(h\), linear convolution is

\[
y[n] = (x*h)[n] = \sum_k x[k]h[n-k].
\]

There is no conjugation in convolution, so the same definition applies to real and complex samples. For finite nonempty inputs of lengths \(N\) and \(M\), the full linear result has length \(N+M-1\).

`convolution::direct` evaluates the definition in \(O(NM)\). `convolution::fft` chooses a power-of-two transform length \(K\ge N+M-1\), zero-pads both operands, computes

\[
Y_k = X_k H_k,
\]

and applies the inverse transform. The padding converts the FFT's circular convolution into linear convolution.

For equal-length sequences, circular convolution is

\[
y[n] = \sum_{k=0}^{N-1}x[k]h[(n-k)\bmod N].
\]

`convolution::circular` evaluates the periodic sum directly. `convolution::circular_bluestein` transforms exactly \(N\) samples with Bluestein, so arbitrary circular lengths do not change period through padding.

### Overlap-add

Split the input into blocks of \(L\) new samples,

\[
x[n] = \sum_b x_b[n-bL].
\]

By linearity,

\[
x*h = \sum_b (x_b*h)[n-bL].
\]

`convolution::overlap_add` FFT-convolves each block with the fixed kernel, then adds overlapping tails at shifted positions. `block_size` is a caller-visible algorithm parameter; the library does not auto-tune it.

### Overlap-save

For kernel length \(M\) and FFT length \(K\ge M\), overlap-save consumes

\[
L = K-M+1
\]

new samples per block. Each transform input contains the previous \(M-1\) samples followed by \(L\) new samples. Circular aliasing is confined to the first \(M-1\) outputs, which are discarded.

`OverlapSavePlan<T>` stores the fixed kernel spectrum and \(M-1\)-sample history. The transform length is explicit and is a power of two because this implementation deliberately uses `Radix2Plan`.

### Streaming direct convolution

`StreamingDirect<T>` keeps the last \(M\) input samples in a circular history buffer and evaluates the defining sum once per new sample. `flush()` feeds the \(M-1\) trailing zeros required to recover the tail of a finite convolution.

### Uniform partitioned convolution

For a long fixed kernel, write

\[
h[n] = \sum_{p=0}^{Q-1} h_p[n-pP]
\]

with partitions of length \(P\). If \(X_b\) is the FFT of input block \(b\) and \(H_p\) the precomputed FFT of kernel partition \(p\), the frequency-domain delay line forms

\[
Y_b[k] = \sum_{p=0}^{Q-1} X_{b-p}[k]H_p[k].
\]

One inverse FFT produces the current block contribution and its second half is overlap-added into the next block. `partition_size` remains explicit because latency, memory, and arithmetic work are properties of the method.

### Convolution numerical behavior

All convolution forms support binary32/binary64 real and complex samples. FFT-based real convolution performs arithmetic in the corresponding complex format and returns the real component after inversion. Blocking and partitioning alter floating-point operation order, so correctness tests use format-specific tolerances.

## Correlation

The library uses

\[
r_{xy}[\ell] = \sum_n x[n]\overline{y[n-\ell]}.
\]

The returned vector stores lag \(\ell\) at index

\[
\ell + (|y|-1),
\]

so zero lag is at index \(|y|-1\).

`correlation::cross_direct` evaluates this sum directly. Reversing and conjugating the second sequence gives

\[
r_{xy} = x * \operatorname{reverse}(\overline{y}),
\]

which is used by `correlation::cross_fft`. `auto_direct` and `auto_fft` set \(y=x\).

Per-lag normalized cross-correlation is

\[
\rho_{xy}[\ell] =
\frac{r_{xy}[\ell]}
{\sqrt{\left(\sum_{n\in I_\ell}|x[n]|^2\right)
       \left(\sum_{n\in I_\ell}|y[n-\ell]|^2\right)}}.
\]

`normalized_direct` and `normalized_fft` use the same overlap-energy denominator. A zero-energy overlap maps to zero. This is energy normalization, not biased/unbiased lag-count normalization.

## Remaining algorithm families

The detailed formulas for the rest of the completed v1 catalog are separated by subject rather than duplicated here:

- [`transforms.md`](transforms.md): DCT-I–IV, DST-I–IV, Hartley, Walsh-Hadamard, MDCT/IMDCT, and window functions;
- [`fir-multirate.md`](fir-multirate.md): FIR execution, linear-phase reductions, windowed-sinc and Parks-McClellan/Remez design, interpolation/decimation, rational conversion, and polyphase forms;
- [`iir.md`](iir.md): Direct Form I/II and transposed realizations, SOS, Butterworth/Chebyshev/elliptic/Bessel prototypes, frequency transformations, and bilinear digital design;
- [`spectral.md`](spectral.md): STFT/ISTFT, periodogram/Welch/CSD/spectrogram density, Goertzel, Hilbert transform, analytic signal, and envelope extraction.

## Numerical conventions

The repository is an algorithm library rather than a benchmark framework. No method silently substitutes another implementation. Stateful objects own only state intrinsic to their algorithm; benchmarking, timing, statistics, campaign orchestration, and reporting remain outside this repository.

Cross-library scalar, empty-input, exception, size-safety, frequency-unit, and compatibility rules are specified in [`contracts.md`](contracts.md).
