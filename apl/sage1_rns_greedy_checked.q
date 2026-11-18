/ SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
/ Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
/ Checked copy of sage1_rns_greedy.q. Change: line 1 made a comment.
/ ===========================================================================
/ SAGE-1 RNS (Residue Number System) Norm & Greedy Neural Network (q / kdb+)
/ Zero-Dependency Pure q Implementation
/ Integrates with sage1_riemannian_jacobian.q & Wolfram Flat CSR Pipe
/ ============================================================================

\P 12

/ --- 1. Core Math & Modular Arithmetic Helpers ---
dot:{sum x*y}
norm2:{sqrt sum x*x}
unitSphere:{[x] n:norm2 x; $[n=0f; x; x%n]}
clip:{[lo;hi;v] lo|hi&v}

/ Extended Euclidean Algorithm: returns (g; x; y) with a*x + b*y = g
egcd:{[a;b]
  $[b=0; (a;1;0);
    [r:.z.s[b; a mod b]; (r 0; r 2; (r 1) - (a div b)*r 2)]]
 }

/ Modular inverse a^{-1} mod m
modInv:{[a;m]
  r:egcd[((a mod m)+m) mod m; m];
  if[1<>r 0; '"Moduli must be pairwise coprime"];
  ((r 1) mod m + m) mod m
 }

/ --- 2. RNS (Residue Number System) Norm ---
/ Pairwise coprime prime moduli (dynamic range prod M ~ 9.4e11)
RNS_MODULI:997 991 983 977j
RNS_SCALE:10000f

/ Encode float vector x into RNS residue matrix (count[M] x count[x])
rnsEncode:{[scale;M;x]
  xq:"j"$x*scale;
  {[xq;m] ((xq mod m)+m) mod m}[xq] each M
 }

/ Chinese Remainder Theorem reconstruction from residue vector r over moduli M
crtReconstruct:{[M;r]
  P:prd M;
  Mi:P div M;
  yi:modInv'[Mi; M];
  (sum r * Mi * yi) mod P
 }

/ Exact RNS Squared Norm: computes sum(x_i^2) in each residue channel,
/ reconstructs exact integer sum via CRT, and rescales to float
rnsNormSq:{[scale;M;x]
  R:rnsEncode[scale;M;x];
  / Carry-free per-modulus inner product: sum(r_k^2) mod m_k
  chanSq:{[r;m] (sum (r*r) mod m) mod m}'[R; M];
  exactInt:crtReconstruct[M; chanSq];
  ("f"$exactInt) % (scale * scale)
 }

/ RNS L2 Norm & Unit-Sphere Projection
rnsNorm:{[scale;M;x] sqrt rnsNormSq[scale;M;x]}

rnsUnitSphere:{[scale;M;x]
  n:rnsNorm[scale;M;x];
  $[n=0f; x; x % n]
 }

/ --- 3. Greedy Neural Network (Geodesic Matching Pursuit + Layer-Wise Fit) ---
/ Geodesic SLERP on S^(n-1)
slerp:{[u;v;lam]
  un:unitSphere u; vn:unitSphere v;
  th:acos clip[-1f;1f;dot[un;vn]];
  $[th<1e-7; unitSphere ((1f-lam)*un)+lam*vn;
    (((sin (1f-lam)*th)%sin th)*un) + ((sin lam*th)%sin th)*vn]
 }

/ Greedy atom selection: picks the dictionary/weight atom maximizing
/ RNS-normalized alignment <residual, atom>
greedySelectAtom:{[scale;M;atoms;target;state]
  res:target - state;
  rNormed:rnsUnitSphere[scale;M;res];
  scores:{[r;a] abs dot[r; unitSphere a]}[rNormed] each atoms;
  bestIdx:first where scores = max scores;
  `idx`atom`score!(bestIdx; atoms bestIdx; scores bestIdx)
 }

/ Single Greedy Neural Layer:
/ 1) Greedily selects best atom from layer codebook
/ 2) Applies ReLU-gated residual step + RNS unit-sphere normalization
/ 3) Evolves state via geodesic SLERP toward greedy target
greedyLayerStep:{[scale;M;atoms;target;lam;state]
  pick:greedySelectAtom[scale;M;atoms;target;state];
  a:rnsUnitSphere[scale;M;pick`atom];
  / Greedy signed projection coefficient
  alpha:dot[target - state; a];
  / Non-linear greedy activation (clipped soft-threshold on tangent)
  act:0f | alpha * a;
  cand:rnsUnitSphere[scale;M; state + act];
  nextState:slerp[state; cand; lam];
  err:norm2[target - nextState];
  `state`atom_idx`alpha`rns_norm`res_err!(
    nextState; pick`idx; alpha; rnsNorm[scale;M;nextState]; err)
 }

/ Multi-Layer Greedy Neural Network Forward Pass
/ Greedily fits & propagates state across L layers without backpropagation
greedyNetForward:{[scale;M;layerCodebooks;target;lam;x0]
  L:count layerCodebooks;
  initTrace:();
  stepFn:{[scale;M;codebooks;target;lam;acc;l]
    st:acc 0; tr:acc 1;
    res:greedyLayerStep[scale;M;codebooks l;target;lam;st];
    row:([] layer:enlist l; atom_idx:enlist res`atom_idx;
            alpha:enlist res`alpha; rns_norm:enlist res`rns_norm;
            res_err:enlist res`res_err);
    (res`state; tr,row)
  }[scale;M;layerCodebooks;target;lam];
  out:stepFn/[(rnsUnitSphere[scale;M;x0]; initTrace); til L];
  `final_state`trace!(out 0; out 1)
 }

/ --- 4. Self-Contained Verification Harness ---
runRnsGreedyTest:{[]
  x:unitSphere 0.42 0.18 0.65 0.31 0.52f;
  rn:rnsNorm[RNS_SCALE; RNS_MODULI; x];
  fn:norm2 x;
  -1 "RNS Norm vs Float Norm2: ",string[rn]," vs ",string[fn]," (diff=",string[abs rn-fn],")";

  / 4-layer greedy neural network with 5D unit-sphere atom codebooks
  cb0:(1 0 0 0 0f; 0 1 0 0 0f; 0 0 1 0 0f; 0 0 0 1 0f; 0 0 0 0 1f);
  cb1:(0.7 0.7 0 0 0f; 0 0.7 0.7 0 0f; 0 0 0.7 0.7 0f; 0 0 0 0.7 0.7f);
  cb2:(0.5 0.5 0.5 0 0f; 0 0.5 0.5 0.5 0f; 0 0 0.5 0.5 0.5f);
  cb3:(0.6 0.2 0.6 0.2 0.4f; 0.2 0.6 0.2 0.6 0.4f);
  codebooks:(cb0; cb1; cb2; cb3);

  target:unitSphere 0.8 0.1 0.5 0.2 0.3f;
  res:greedyNetForward[RNS_SCALE; RNS_MODULI; codebooks; target; 0.45f; x];
  show res`trace;
  res
 }

/ Execute self-test on load
runRnsGreedyTest[]
