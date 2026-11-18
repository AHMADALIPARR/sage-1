# SAGE-1 — native CSR grammar / spherical array kernel

Standalone C11 module for inclusion in `uor-r4`. It does **not** modify the existing Rust qualification path (`crates/uor-r4-core/src/native_support.rs`). The repository's `src/lib.rs` is not its compiled crate entry point; use a separately gated FFI boundary after review.

Build and test:

```sh
cc -std=c11 -O2 -Wall -Wextra -Werror -pedantic -fsanitize=address,undefined -fno-omit-frame-pointer native/sage1/sage1.c native/sage1/test_sage1.c -lm -o /tmp/uor-sage1-test
/tmp/uor-sage1-test
```

Algorithm:

1. Sparse CSR counts determine conditional row distributions (no learned neural weight tensor).
2. `q = P^T (x*x)` (O(edges + states)); `y = sqrt(q)`.
3. `out = SLERP(x,y,lambda)` preserves unit-sphere geometry.
4. A grammar engine supplies `legal` masks; `sage1_select` rejects dead ends.
5. Call `sage1_condition` *after updating the parser state* to retain the selected branch.

`lambda=1` reproduces Markov propagation exactly after squaring; intermediate steps differ. Counts supplied by a grammar learner and stack updates from Prolog are not yet connected. A static CSR model is not itself a full pushdown grammar or a conversational model. Float64 is a research-side numerical boundary, not the existing integer-artifact serving pathway.

No external dependencies beyond C11 and libm. Arrays are caller-owned. Inputs must not alias outputs. Rows must be nonempty, entries must have positive counts, row totals must fit uint64, and every state must have at least one successor. The decoder's deterministic argmax resolves ties by smallest state ID.
