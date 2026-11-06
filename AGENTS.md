# Signal Processing Research Task

This repository contains a method-explicit, header-only C++23 signal processing library. The goal is correctness, numerical reliability, performance measurement, and reproducible research-quality documentation—not merely passing one build.

## Working branch
Modify only `ai/localci-dsp-research`. Never merge or push to `main` or remove an algorithm to make tests pass. Use `EddyTodd/local-ci-coordinator` for all benchmarking and machine measurement infrastructure. Keep build outputs and benchmarks out of this tracked source checkout.

## Correctness gates
- Reference independent mathematical definitions and cross-method equivalence for float32 and float64.
- Evaluate prime/composite/radix-two lengths, non-power-of-two transforms, short and empty/invalid inputs, inverse normalization, filter states, boundary indexing, NaN/Inf and overflow contracts where applicable.
- Exercise metamorphic properties, conservation laws, time-frequency reconstruction, stable IIR poles and correct filter gain.
- Verify assertions remain active in all build modes. Never substitute a quick smoke build for tests.
- Run Release CTest and ASan/UBSan locally, with test counts and exact compiler/build flags in artifacts.
- A PASS on finite tests is not a numerical proof, stability guarantee for all parameters, or certification of optimality.

## Performance gates
- Benchmark on the same machine with pinned source, compiler version, ISA, CPU model, fixed seeds, repeated timed samples and spread.
- Compare only matched precision, method setup, allocation behavior, threading and input sizes. Publish complete raw samples and the median and variability, not cherry-picked minimum.
- Include DFT, iterative/recursive radix variants, prime and mixed-length transforms, preplanned kernels, convolution, FIR/DCT, and IIR/STFT as available.
- Cross-compiler, float32, memory/throughput, adversarial numerical inputs and third-party baseline comparisons are further milestones. Do not claim them done without evidence.

## Website
- Publish an understandable project overview, formulas and sampling conventions, comparison tables, graphics and interactive charts.
- Only ingest real benchmark data published by Local CI; never fabricate measurements. Clearly mark pending/illustrative figures.
- Label compiler and CPU, size, precision, timing regime, uncertainty, source SHA, and comparison limitations. Preserve accessibility and mobile usability.

## Automation contract
The Local CI task ID is `signal-processing` and the project ID is `signal_processing`. Read the exact report before each iteration, fix evidenced failures, and update `.localci/state.json`. When all acceptance gates are demonstrably satisfied, set state COMPLETE and update coordinator task status only after CI evidence. If not verified, keep ACTIVE. A final handoff-only commit changing `automation/handoff.json` must be the last change of a continuation cycle; no other file belongs in it. Do not create the handoff before completing source changes and report review.
