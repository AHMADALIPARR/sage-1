<!-- SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0 -->
<!-- Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab -->

# SAGE-1

[![CI](https://github.com/AHMADALIPARR/sage-1/actions/workflows/ci.yml/badge.svg)](https://github.com/AHMADALIPARR/sage-1/actions/workflows/ci.yml)
[![License: AGPL-3.0 OR SnapKitty Opaque 1.0](https://img.shields.io/badge/license-AGPL--3.0--only%20OR%20SnapKitty%20Opaque%201.0-blue)](LICENSE)
![C11](https://img.shields.io/badge/C-C11-555555?logo=c)
![Rust](https://img.shields.io/badge/Rust-runtime-b7410e?logo=rust)
![Prolog](https://img.shields.io/badge/SWI--Prolog-grammar-e61b23)
![Futhark](https://img.shields.io/badge/Futhark-array%20backend-5f021f)
![TinyAPL](https://img.shields.io/badge/TinyAPL-array%20engine-2a6db0)
![q](https://img.shields.io/badge/q%2Fkdb%2B-attention%20head-1e78b3)

**Sparse grammar transitions evolved on the unit sphere.** SAGE-1 replaces
attention weights with integer transition counts: a sparse CSR count matrix
defines a row-stochastic operator `P`, the state lives on the unit sphere as
`x` with probabilities `x ⊙ x`, and every step moves along a great circle:

```
x_{t+1} = SLERP( x_t , sqrt( Pᵀ (x_t ⊙ x_t) ) , λ )
```

There are no QKV projections and no neural weight tensors. A grammar decides
which transitions are legal, decoding is a deterministic argmax over the legal
set, and the state is collapsed onto the chosen branch after the parser has
moved.

## Features

- **C11 kernel** (`native/sage1`): CSR validation, count-normalized propagation
  in `O(edges + states)`, geodesic SLERP that stays in the nonnegative orthant
  of the sphere, masked deterministic decoding with dead-end detection, and
  state conditioning. No dependencies beyond libm; all buffers caller-owned.
- **Rust runtime** (`rust/sage1`): safe bindings with every length checked
  before the FFI call, an allocation-free `step_into`, and a generic
  grammar-constrained decoder that backtracks over dead ends.
- **Prolog grammar** (`prolog/`): a SWI-Prolog foreign library over the kernel
  plus a pushdown grammar in Prolog that produces the legal masks and the stack
  updates; dead ends become Prolog backtracking.
- **Futhark array backend** (`futhark/`): data-parallel propagation (segmented
  CSR scatter-reduce), sqrt + SLERP step, fused multi-step evolution, decoding
  and conditioning, built as a C library for the `c` and `multicore` backends.
- **Step variants** (`native/sage1/variants`): a temperature-scaled
  log-ratio step that replaces the SLERP angle, and a low-rank residual
  projection between CSR propagation and the grammar mask, fitted
  deterministically from training transitions.
- **Perplexity benchmark** (`bench/`): masked Markov softmax, the SLERP step
  and both variants scored on a held-out split under the same grammar mask,
  with the C grammar port checked against the Prolog grammar.
- **Array engine** (`apl/`): the TinyAPL engine with columnar CSR and Key
  scatter-reduce, a q/kdb+ Riemannian attention head with step Jacobians, and an
  RNS exact-norm greedy layer network, joined by a connector script.

## Quick start

```sh
cmake -S . -B build -DSAGE1_SANITIZE=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Everything optional is a CMake switch:

```sh
cmake -S . -B build -DSAGE1_PROLOG=ON -DSAGE1_FUTHARK=ON -DSAGE1_TINYAPL=ON -DSAGE1_BENCH=ON
```

Rust:

```sh
cd rust/sage1 && cargo test
```

The original sanitizer build from `native/sage1/README.md` also works as-is:

```sh
cc -std=c11 -O2 -Wall -Wextra -Werror -pedantic -fsanitize=address,undefined -fno-omit-frame-pointer native/sage1/sage1.c native/sage1/test_sage1.c -lm -o /tmp/uor-sage1-test
/tmp/uor-sage1-test
```

## The algorithm

1. Row `i` of the CSR matrix holds counts `i → j`; `P_ij = count_ij / Σ_k count_ik`.
2. `p = x ⊙ x` is a probability vector because `‖x‖₂ = 1`.
3. `q = Pᵀ p`, then `y = sqrt(q)` is again a unit vector.
4. `out = SLERP(x, y, λ)` moves a fraction `λ` of the way along the great
   circle. `λ = 1` reproduces the Markov chain exactly after squaring; `λ = 0`
   stays put; intermediate values keep `‖out‖₂ = 1`.
5. The grammar supplies a `legal` mask. `sage1_select` returns the legal state
   with the largest `x_i²` (ties → smallest index) and its probability
   renormalized over the legal set, or `SAGE1_DEAD_END` when no legal state has
   mass.
6. The parser state is updated with the selected token, and only then
   `sage1_condition` collapses `x` onto it.

## C API

```c
#include "sage1.h"

typedef struct {
    size_t states, edges;
    const size_t   *row_ptr;   /* states + 1 */
    const uint32_t *col_idx;   /* edges */
    const uint64_t *counts;    /* edges */
} sage1_csr;

sage1_status sage1_validate (const sage1_csr *m);
sage1_status sage1_propagate(const sage1_csr *m, const double *p, double *q);
sage1_status sage1_geodesic (const double *x, const double *y, size_t n, double lambda, double *out);
sage1_status sage1_step     (const sage1_csr *m, const double *x, double lambda,
                             double *out, double *q, double *y);
sage1_status sage1_select   (const double *x, const uint8_t *legal, size_t n,
                             size_t *selected, double *probability);
sage1_status sage1_condition(size_t selected, size_t n, double *x);
```

Status codes: `SAGE1_OK`, `SAGE1_BAD_INPUT`, `SAGE1_BAD_MATRIX`,
`SAGE1_DEAD_END`, `SAGE1_NUMERIC`. Rows must be nonempty with positive counts,
vectors must be finite and unit-norm in the nonnegative orthant, and inputs
must not alias outputs.

## Step variants

```c
#include "sage1_variants.h"

/* A: q = Pᵀ(x ⊙ x); out = q^(1/(2τ)) / ‖q^(1/(2τ))‖ */
sage1_status sage1_logratio_step(const sage1_csr *m, const double *x, double tau,
                                 double *out, double *q);
/* B: y = sqrt(Pᵀ(x ⊙ x)); z = y + U(Vᵀy); out = SLERP(x, |z|/‖z‖, λ) */
sage1_status sage1_lowrank_step(const sage1_csr *m, const double *x, double lambda,
                                const double *U, const double *V, size_t rank,
                                double *out, double *y, double *z, double *w);
sage1_status sage1_lowrank_fit(const double *S, const double *T, size_t n,
                               size_t rank, double *U, double *V);
```

**A — temperature-scaled log-ratio step.** The fixed angle `λ` is replaced by
a step that scales every log-ratio of the propagated target by `1/τ`:
`log(out_i² / out_j²) = (1/τ)·log(q_i / q_j)`. Zero entries stay zero, the
output is a nonnegative unit vector, and `τ = 1` is `sage1_step` at `λ = 1`.
It is computed with a max shift in the log domain, so small `τ` does not
overflow.

**B — low-rank residual projection.** After CSR propagation the amplitude
`y` gets a rank-`r` correction `U(Vᵀy)`, is folded back into the nonnegative
orthant with `|·|` (probabilities `z_i²` are unchanged) and renormalized, and
the SLERP step and grammar mask follow as before. `sage1_lowrank_fit` takes
the propagated amplitude rows `S` (row `k` = `sqrt(P[k,:])`) and target rows
`T`, solves `S X = T − S` by Gaussian elimination, and keeps the top-`r` left
singular vectors `U_r` of `(T − S)ᵀ` (cyclic Jacobi, fixed sweep order, sign
fixed so each vector's largest entry is positive): `U = U_r`, `V = X U_r`.
That is the rank-`r` reduced-rank regression of the residual on `S`; rank
`n` reproduces every target row exactly. `rank = 0` is `sage1_step`.

## Rust

```rust
use sage1::{Csr, select, condition};
use sage1::grammar::brackets::{self, Brackets};

let m = Csr::from_triples(3, &[(0, 1, 3), (0, 2, 1), (1, 2, 1), (2, 0, 1)])?;
let x = m.step(&[1.0, 0.0, 0.0], 0.5)?;          // SLERP step
let s = select(&x, &[false, true, true])?;       // masked decode
let x = condition(s.index, m.states())?;         // collapse

// Grammar-constrained decoding with backtracking over dead ends
let d = sage1::decode(&brackets::csr()?, &Brackets::new(4), brackets::BOS, 1.0, 16)?;
```

`build.rs` compiles `native/sage1/sage1.c` with the system C compiler (`$CC`,
`$AR`, extra flags via `$SAGE1_CFLAGS`); the crate has no dependencies.
Implement `Grammar` (`legal`, `advance`, `is_final`) for any parser to drive
`decode`.

## Prolog

```prolog
?- use_module(prolog/sage1), use_module(prolog/sage1_grammar).
?- demo_csr(Csr), initial_parser(4, P0),
   decode(Csr, P0, 1.0, 16, Tokens, Trace).
Tokens = [3, 1, 5, 2, 4, 6]      % [ ( a ) ] </s>
```

`sage1_pl.c` exposes `sage1_validate/2`, `sage1_propagate/3`, `sage1_step/4`,
`sage1_geodesic/4`, `sage1_select/3` (`selected(I, P)` or `dead_end`) and
`sage1_condition/3`, with matrices as `csr(States, RowPtr, ColIdx, Counts)`.
Kernel errors raise `error(sage1_error(Status), _)`. Point `SAGE1_PL_LIB` at
the built `sage1_pl.so`, or put it on the `foreign/1` search path.

The demo grammar is a pushdown automaton over `<s> ( ) [ ] a </s>`: brackets
must match the stack top, atoms cannot repeat, depth is bounded and `</s>` is
legal only on an empty stack. With the demo counts the greedy path `[ a`
dead-ends (after `a` the corpus only continues with `)` or `</s>`, both
illegal inside `[`); `sage1_select` reports it, Prolog backtracks and the
decoder emits `[ ( a ) ] </s>`. A separate DCG for the same language checks
every decoded sentence. The Rust decoder produces the same tokens, the same
probabilities `0.5, 1, 0.8, 1, 1, 2/3` and the same dead end.

## Futhark

`futhark/sage1.fut` provides the entries `propagate`, `step`, `evolve`
(fused multi-step loop), `geodesic`, `select` and `condition`, each returning
the `sage1_status` code first. CMake runs `futhark <backend> --library` for
every backend in `SAGE1_FUTHARK_BACKENDS` (default `c;multicore`) and links
`futhark/test_futhark_vs_c.c`, which runs both implementations on the same
pseudo-random CSR matrices up to 20,000 states (propagation, steps at
`λ ∈ {0, 0.25, 0.5, 0.9, 1}`, 20-step evolution, geodesics, masked decoding
and every error status) and requires agreement within `1e-12`.

## Array engine

The `apl/` directory is the array-language layer. Original files are kept
byte-for-byte; the `_checked` copies carry the fixes needed to run them.

| File | Role |
|---|---|
| `sage1_engine.tinyapl` | TinyAPL engine (original), followed by the Wolfram Language `SAGE1GeometryPipe` package that triangulates geometry and exports the flat CSR pipe |
| `sage1_engine_checked.tinyapl` | Runnable TinyAPL 0.12 copy of the engine |
| `sage1_attention_head.q` | q attention head: Riemannian manifold attention and step Jacobians (original) |
| `sage1_rns_greedy.q` | q RNS exact norm and greedy layer network (original) |
| `sage1_attention_head_checked.q`, `sage1_rns_greedy_checked.q` | Checked q copies loaded by the connector |
| `sage1_rns_attention.q` | Connector: CSV pipe → RNS attention → SAGE-1 step → greedy network |

**TinyAPL engine.** `Norm2`, `UnitSphere` and `ProbDecode` are the sphere
primitives, `Slerp` interpolates with a small-angle fallback, and the CSR is a
flat columnar record `⟨n ⋄ row ⋄ col ⋄ w⟩`. `SparsePtMul` computes `Pᵀ p`
without a dense matrix: per-edge flow `p[row] × w`, grouped by target with
Key (`⌸`) and scattered back. `SAGE1Step` composes these into the SAGE-1 step
and checks `‖x‖₂ = 1`, `StateCondition` boosts peak states and re-projects, and
`DetSelect` is the argmax decoder. On the 5-vertex triangulated patch it
reproduces the C kernel: CI runs the checked copy under TinyAPL and compares
`x₁`, `x₂`, both norms and the selected state with `sage1_step`/`sage1_select`.

**Attention head (q).** Weight-free attention on `Sⁿ⁻¹`: weights are
`softmax(−d²/τ)` of squared geodesic distances between the block's states, and
each state moves by the weighted mean of its log-map tangents, mapped back
with the exponential map (a one-step Fréchet mean update). Each block then
applies the SAGE-1 SLERP step to every state, composes the step Jacobian
`(1−λ)I + λ·J_proj(√q)·diag(1/2√q)·Pᵀ·diag(2x)` into an accumulated Jacobian,
and records `depth`, `mean_norm` and `jac_frob`; `fwdBlockRec` recurses over
the blocks. `loadWolframCSV` reads `sage1_csr.csv`, `sage1_peaks.csv` and
`sage1_init.csv` from the Wolfram pipe and converts its 1-based indices to
0-based.

**RNS greedy network (q).** `rnsNormSq` quantizes a vector at scale `10⁴`,
squares and sums it independently in four residue channels (moduli 997, 991,
983, 977, product ≈ 9.4·10¹¹) and rebuilds the exact integer with the Chinese
Remainder Theorem, so the norm of the quantized vector carries no
floating-point rounding. `rnsUnitSphere` projects with that norm. The greedy
network needs no backpropagation: each layer picks the codebook atom best
aligned with the RNS-normalized residual `target − state`, takes a ReLU-gated
step along it, re-projects, and SLERPs toward that candidate; the trace
records the atom, the coefficient, the RNS norm and the residual error per
layer.

**Connector data flow** (`q apl/sage1_rns_attention.q [pipe_dir]`, run from the
repository root):

```
Wolfram flat CSV pipe ── loadWolframCSV ──▶ (n, CSR, peaks, x₀)      or the built-in 5-vertex patch
        │
        ▼  per block (λ, τ) = (0.30, 0.50), (0.35, 0.40), (0.40, 0.30)
RNS-projected Riemannian attention  ── unitSphere bound to rnsUnitSphere
        │
        ▼
SAGE-1 SLERP step per state + step Jacobian  ── float unitSphere (same step as C / TinyAPL)
        │
        ▼
final attention read-out ──▶ greedy layer network toward its SAGE-1 target
                             (codebook: rows of P, square-rooted onto the sphere)
        │
        ▼
combined trace: one row per attention block, one row per greedy layer
```

Both checked scripts define `dot`, `norm2`, `unitSphere`, `clip` and `slerp`
with identical bodies, so they load in sequence. Because q resolves globals at
call time, the connector rebinds `unitSphere` to the RNS projection only while
attention runs (`.sra.with`) and restores the float projection for the SAGE-1
step, whose unit-norm check is `1e-8`.

## Benchmarks

```sh
cmake -S . -B build -DSAGE1_BENCH=ON && cmake --build build
./build/sage1_bench            # --sentences N --seed S --max-depth D --out FILE
```

**Corpus.** 20,000 sentences sampled (splitmix64, seed 20261118) from a
stochastic version of the Prolog grammar's language with stack depth ≤ 3:
each nesting level stops or adds an atom or a `( … )` / `[ … ]` group, with
probabilities that depend on the depth and the enclosing bracket. The first
16,000 sentences (185,734 tokens) are training data and the last 4,000
(45,237 tokens) are held out. `bench/check_grammar.pl` checks the C grammar
port against `legal_mask/2` on all 105 parser states and parses 2,000
generated sentences with `sentence//0`.

**Protocol.** Every model gets the same add-one bigram counts from the
training split as its CSR matrix. At each held-out position the state is
`e_k` for the previous token `k` (as after `sage1_condition`), the model
produces amplitudes, and the probability of token `j` is `x_j²` over the
legal set — what `sage1_select` scores. The softmax baseline is
`p_j ∝ exp(log P[k][j] / T)` over the legal set. Variant B's targets are the
masked choice-model fit on the training split: for each `k`, the weights `w`
that maximize `Σ log(w_next / Σ_legal w)` with one add-one pseudo-event per
pair, found with 2,000 MM iterations; `T[k] = sqrt(w / Σw)`. Margin is the
mean of top-1 minus top-2 probability over positions with at least two legal
tokens; a mask violation is a position where `sage1_select` returns an
illegal token or the held-out token is illegal.

|---|---|---:|---:|---:|---:|
| Generator (true distribution) | — | 3.0809 | — | — | — |
| Masked Markov softmax | T = 0.50 | 3.8655 | 0.3557 | 39.16% | 0 |
| Masked Markov softmax | T = 0.75 | 3.4739 | 0.2400 | 39.16% | 0 |
| Masked Markov softmax | T = 1.00 | 3.3701 | 0.1763 | 39.16% | 0 |
| Masked Markov softmax | T = 1.25 | 3.3399 | 0.1377 | 39.16% | 0 |
| Masked Markov softmax | T = 1.50 | 3.3342 | 0.1124 | 39.16% | 0 |
| Masked Markov softmax | T = 2.00 | 3.3440 | 0.0816 | 39.16% | 0 |
| Masked Markov softmax | T = 2.50 | 3.3594 | 0.0638 | 39.16% | 0 |
| SAGE-1 SLERP step | λ = 0.25 | 5.8970 | 0.4323 | 31.17% | 0 |
| SAGE-1 SLERP step | λ = 0.30 | 5.3032 | 0.4099 | 31.17% | 0 |
| SAGE-1 SLERP step | λ = 0.50 | 4.0278 | 0.2899 | 31.17% | 0 |
| SAGE-1 SLERP step | λ = 0.90 | 3.3347 | 0.1466 | 42.12% | 0 |
| SAGE-1 SLERP step | λ = 1.00 | 3.3701 | 0.1763 | 39.16% | 0 |
| A: log-ratio step | τ = 0.50 | 3.8655 | 0.3557 | 39.16% | 0 |
| A: log-ratio step | τ = 0.75 | 3.4739 | 0.2400 | 39.16% | 0 |
| A: log-ratio step | τ = 1.00 | 3.3701 | 0.1763 | 39.16% | 0 |
| A: log-ratio step | τ = 1.25 | 3.3399 | 0.1377 | 39.16% | 0 |
| A: log-ratio step | τ = 1.50 | 3.3342 | 0.1124 | 39.16% | 0 |
| A: log-ratio step | τ = 2.00 | 3.3440 | 0.0816 | 39.16% | 0 |
| A: log-ratio step | τ = 2.50 | 3.3594 | 0.0638 | 39.16% | 0 |
| B: low-rank residual | rank 0, λ = 1 | 3.3701 | 0.1763 | 39.16% | 0 |
| B: low-rank residual | rank 1, λ = 1 | 3.3774 | 0.3317 | 46.73% | 0 |
| B: low-rank residual | rank 2, λ = 1 | 3.4122 | 0.2471 | 46.65% | 0 |
| B: low-rank residual | rank 3, λ = 1 | 3.2599 | 0.2455 | 45.45% | 0 |
| B: low-rank residual | rank 4, λ = 1 | 3.1519 | 0.2085 | 47.76% | 0 |
| B: low-rank residual | rank 5, λ = 1 | 3.1239 | 0.2010 | 48.06% | 0 |
| B: low-rank residual | rank 6, λ = 1 | 3.1239 | 0.2007 | 48.06% | 0 |
| B: low-rank residual | rank 7, λ = 1 | 3.1239 | 0.2007 | 48.06% | 0 |

**Where the gap closes.** The best softmax temperature in the sweep is
`T = 1.5` at 3.3342. The SLERP step is identical to softmax `T = 1` at
`λ = 1` (3.3701) and comes within 0.0005 of the best softmax at `λ = 0.9`
(3.3347); at `λ ≤ 0.5` the extra mass on the current token raises margin but
costs perplexity (5.8970 at `λ = 0.25`). Variant A matches softmax at
`T = τ` on every row: with a one-hot state it computes the same masked
`q^(1/τ)` distribution, so it closes the gap completely and inherits the
same perplexity/margin trade. Variant B is the only one that goes past
softmax: rank 1 keeps perplexity near softmax `T = 1` (3.3774) with nearly
the margin of `T = 0.5` (0.3317 vs 0.3557), rank 3 (3.2599) beats every
softmax temperature, and ranks 5–7 reach 3.1239 with 48.06% top-1 accuracy,
0.043 above the generator's own 3.0809. Its gain comes from fitting the
residual to the masked likelihood, which plain bigram counts do not see.
No method produced a mask violation.

## Layout

```
native/sage1/      C11 kernel, original test and README (unmodified)
native/sage1/variants/  step variants A (log-ratio) and B (low-rank residual)
tests/c/           extended C tests, README sanitizer command
rust/sage1/        Rust crate: build.rs, FFI, safe API, grammar decoder, tests
prolog/            foreign library, Prolog binding, pushdown grammar, plunit tests
futhark/           Futhark kernels and the C cross-check driver
apl/               TinyAPL engine, q attention head, RNS network, connector
tests/apl/         TinyAPL-vs-C check
tests/variants/    variant A/B checks
bench/             perplexity benchmark, grammar port, Prolog cross-check
cmake/             SWI-Prolog discovery
```

## Build requirements

| Part | Needs |
|---|---|
| C kernel and tests | C11 compiler, CMake ≥ 3.16; ASan/UBSan for `SAGE1_SANITIZE` |
| Rust | Rust ≥ 1.70 and a C compiler |
| Benchmark | C kernel requirements; SWI-Prolog for the grammar cross-check |
| Prolog | SWI-Prolog with headers (`swipl --dump-runtime-variables` locates them) |
| Futhark | `futhark` ≥ 0.25 (CI uses 0.29.1); pthreads for `multicore` |
| TinyAPL | `tinyapl` 0.12 CLI built with GHC 9.10 |
| q scripts | a q/kdb+ interpreter |

## License

Dual licensed: GNU Affero General Public License v3.0 only
([LICENSE-AGPL-3.0](LICENSE-AGPL-3.0)) **or** SnapKitty Opaque License 1.0
([LICENSE-SnapKitty-Opaque-1.0](LICENSE-SnapKitty-Opaque-1.0)), at your
option. See [LICENSE](LICENSE).

`SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0`

Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
