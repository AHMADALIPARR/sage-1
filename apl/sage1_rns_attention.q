/ SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
/ Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
/ ============================================================================
/ SAGE-1 RNS attention connector (q / kdb+)
/ Wolfram flat CSV pipe -> RNS-projected Riemannian attention
/   -> SAGE-1 SLERP step + step Jacobian (float, matches the C kernel)
/   -> greedy layer network on the attention output.
/
/ Run from the repository root:
/   q apl/sage1_rns_attention.q              / built-in 5-vertex fixture
/   q apl/sage1_rns_attention.q <pipe_dir>   / sage1_csr.csv, sage1_peaks.csv, sage1_init.csv
/
/ Name clashes: both checked scripts define dot, norm2, unitSphere, clip and
/ slerp with identical bodies, so loading them in sequence is order-safe.
/ Every function resolves `unitSphere` globally at call time, so the connector
/ swaps the projection only for the attention stage (.sra.with) and restores
/ the float projection for the SAGE-1 step, whose unit-norm check is 1e-8.
/ Neither original is modified; both checked copies run their own self-tests
/ when loaded.
/ ============================================================================

\l apl/sage1_attention_head_checked.q
\l apl/sage1_rns_greedy_checked.q

.sra.floatUnit:unitSphere;
.sra.rnsUnit:{[x] rnsUnitSphere[RNS_SCALE; RNS_MODULI; x]};

/ Call f . args with `unitSphere` bound to proj, restoring it afterwards.
.sra.with:{[proj;f;args]
  old:get `unitSphere;
  `unitSphere set proj;
  r:.[f; args; {[old;e] `unitSphere set old; 'e}[old]];
  `unitSphere set old;
  r
 }

/ Block recursion: RNS attention, then the float SAGE-1 step per state,
/ Jacobian composition, one trace row per block.
.sra.fwd:{[blocks;X;jAcc;depth;trace]
  if[depth >= count blocks; :(X; jAcc; trace)];
  blk:blocks depth;
  n:blk`n; csr:blk`csr; lam:blk`lam; tau:blk`tau;
  Xattn:.sra.with[.sra.rnsUnit; riemannAttn; (tau; X)];
  Xnext:{[n;csr;lam;x]
    xt:unitSphere x;
    tgt:unitSphere sqrt 0f|sparsePtMul[n;csr;xt*xt];
    xn:slerp[xt; tgt; lam];
    if[1e-8 < abs 1f - norm2 xn; '"Unit-sphere invariant violated"];
    xn
  }[n;csr;lam] each Xattn;
  jNext:sage1Jacobian[n; csr; Xattn 0; lam] mmu jAcc;
  rec:([] stage:enlist `attn_step; idx:enlist depth;
          attn_rns_norm:enlist avg rnsNorm[RNS_SCALE; RNS_MODULI] each Xattn;
          attn_float_norm:enlist avg norm2 each Xattn;
          mean_norm:enlist avg norm2 each Xnext;
          jac_frob:enlist sqrt sum over jNext*jNext);
  .z.s[blocks; Xnext; jNext; depth+1; trace,rec]
 }

/ Transition atoms: row i of P (= P^T e_i), square-rooted onto the sphere.
.sra.atoms:{[n;csr] {[n;csr;i] unitSphere sqrt 0f|sparsePtMul[n;csr;"f"$i=til n]}[n;csr] each til n};

.sra.run:{[pipe]
  n:pipe`n; csr:pipe`csr; pks:pipe`peaks; x0:pipe`x0;
  / second state: x0 with the peaks boosted, back on the sphere
  X0:(x0; unitSphere x0 * 1f + 0.25f * (til n) in pks);
  blocks:(
    `n`csr`lam`tau!(n; csr; 0.30f; 0.50f);
    `n`csr`lam`tau!(n; csr; 0.35f; 0.40f);
    `n`csr`lam`tau!(n; csr; 0.40f; 0.30f));
  res:.sra.fwd[blocks; X0; "f"$eye n; 0; ()];
  / final attention read-out, then greedy fit toward its SAGE-1 target
  Xa:.sra.with[.sra.rnsUnit; riemannAttn; (0.30f; res 0)];
  xa:Xa 0;
  target:unitSphere sqrt 0f|sparsePtMul[n; csr; xa*xa];
  atoms:.sra.atoms[n; csr];
  g:greedyNetForward[RNS_SCALE; RNS_MODULI; (count blocks)#enlist atoms; target; 0.45f; xa];
  gt:update stage:`greedy from select idx:layer, atom_idx, alpha, rns_norm, res_err from g`trace;
  combined:(res 2) uj gt;
  -1 "SAGE-1 RNS attention: combined trace";
  show combined;
  -1 "final greedy state: ",-3!g`final_state;
  `blocks`jacobian`greedy`trace!(res 0; res 1; g; combined)
 }

/ Built-in fixture: the 5-vertex patch from both self-tests (0-based CSR).
.sra.fixture:{[]
  n:5;
  csr:([] row:0 0 1 1 1 2 2 3 3 4 4i;
          col:0 1 0 1 2 1 3 2 4 3 4i;
          w:0.5 0.5 0.3 0.4 0.3 0.5 0.5 0.5 0.5 0.5 0.5f);
  `n`csr`peaks`x0!(n; csr; 0 2 4; unitSphere n#1f)
 }

.sra.result:.sra.run $[count .z.x; loadWolframCSV first .z.x; .sra.fixture[]];
