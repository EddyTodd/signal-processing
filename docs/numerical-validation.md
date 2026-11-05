# Numerical validation and research-grade evidence

This library is actively evaluated; passing included tests does not prove global error bounds or superiority to external DSP libraries.

## Validation tiers

1. Independent textbook O(N²) binary64 DFT and linear-convolution oracles.
2. Parseval's theorem, circular-shift phase, inverse round-trip and DCT-II equivalence.
3. Existing algorithm and filter design regression suites.
4. Independent AddressSanitizer and UndefinedBehaviorSanitizer build where available.
5. External Local CI microbenchmarks with compiler, flags, CPU, ISA, samples, spread, repetitions and allocation regime documented.

## Crucial Release-build contract

The test binaries use C assert statements. A normal Release configuration defines NDEBUG, disabling those assertions. CMake now explicitly undefines NDEBUG on test binaries for every build type, and the numerical test deliberately fails compilation if this policy is violated.

## Limitations

Textbook binary64 oracles are algorithmically independent, but not exact arbitrary-precision reference calculations. Finite seeded tests are regression evidence, not formal proof. Performance data must be measured and attributed to a particular machine and build. Never present illustrative numbers as measured results.
