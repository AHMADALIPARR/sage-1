/ SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
/ Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
/ Checked copy of sage1_attention_head.q. Changes: line 1 made a comment;
/ float identity for mmu (jacProj fallback, self-test J_acc); CSV pipe column
/ `weight` renamed to `w` on load to match sparsePtMul.
/ ===========================================================================
/ SAGE-1 Riemannian Manifold Attention & Jacobian Block Recursive Forward Pass
/ Zero-Dependency q / kdb+ Backend Script
/ Connects to Wolfram Flat CSR Pipe (sage1_geometry_pipe.wl) & TinyAPL Engine
/ Invariant: x in S^(n-1), ||x||_2 = 1
/ ============================================================================

\P 12

/ --- 1. Core Vector & Manifold Primitives on S^(n-1) ---
dot:{sum x*y}
norm2:{sqrt sum x*x}
unitSphere:{[x] n:norm2 x; $[n=0f; x; x%n]}
clip:{[lo;hi;v] lo|hi&v}
outer:{[a;b] a*/:\:b}
eye:{[n] {x=/:x} til n}

/ Sparse CSR transpose-multiply: y = P^T * p using q table grouping
/ csr is a table ([] row:`int$(); col:`int$(); w:`float$()), 0-indexed
sparsePtMul:{[n;csr;p]
  flow:(p csr`row)*csr`w;
  agg:select s:sum flow by col from update flow:flow from csr;
  out:n#0f;
  @[out; key[agg]`col; +; exec s from agg]
 }

/ --- 2. Riemannian Manifold Geometry (LogMap, ExpMap, SLERP) ---
/ Logarithmic map on S^(n-1): Log_x(y) in T_x S^(n-1)
logMap:{[x;y]
  u:unitSphere x; v:unitSphere y;
  c:clip[-1f;1f;dot[u;v]];
  th:acos c;
  $[th<1e-7; v - c*u; (th % sin th) * (v - c*u)]
 }

/ Exponential map on S^(n-1): Exp_x(v) -> S^(n-1)
expMap:{[x;v]
  u:unitSphere x;
  nv:norm2 v;
  $[nv<1e-7; unitSphere u+v; ((cos nv)*u) + ((sin nv)%nv)*v]
 }

/ Geodesic SLERP on S^(n-1)
slerp:{[u;v;lam]
  un:unitSphere u; vn:unitSphere v;
  c:clip[-1f;1f;dot[un;vn]];
  th:acos c;
  $[th<1e-7; unitSphere ((1f-lam)*un)+lam*vn;
    (((sin (1f-lam)*th)%sin th)*un) + ((sin lam*th)%sin th)*vn]
 }

/ --- 3. Geometric Riemannian Manifold Attention ---
/ Pairwise geodesic distance on S^(n-1)
geoDist:{[x;y] acos clip[-1f;1f;dot[unitSphere x; unitSphere y]]}

/ Softmax over a float vector (numerically stabilized)
softmax:{[v] e:exp v - max v; e % sum e}

/ Manifold Attention over a block of unit-sphere state vectors X (k x n)
/ Aggregates in tangent space T_{x_i} S^(n-1) and maps back via ExpMap
riemannAttn:{[tau;X]
  k:count X;
  / Pairwise squared geodesic distance matrix D^2
  d2:{[X;i] {x*x} each geoDist[X i] each X}[X] each til k;
  W:softmax each neg d2 % tau;
  / Tangent-space weighted Fréchet update per token/node i
  {[X;W;i]
    xi:X i; wi:W i;
    tangents:logMap[xi] each X;
    vbar:sum wi * tangents;
    expMap[xi; vbar]
  }[X;W] each til k
 }

/ --- 4. Exact Jacobian Blocks for SAGE-1 & Manifold Projection ---
/ Projection Jacobian J_proj(z) = (I - zhat * zhat^T) / ||z||_2
jacProj:{[z]
  n:count z; nz:norm2 z;
  $[nz<1e-12; "f"$eye n; (eye[n] - outer[z%nz; z%nz]) % nz]
 }

/ Dense materialization of P^T (for small/medium block Jacobian evaluation)
csrToDenseT:{[n;csr]
  M:(n;n)#0f;
  {[M;r] .[M;(r`col;r`row);+;r`w]}/[M;csr]
 }

/ Single-step SAGE-1 Jacobian J_step = d(x_{t+1}) / d(x_t)
/ x_{t+1} = SLERP(x_t, tgt(x_t), lam), tgt = unitSphere(sqrt(P^T (x_t * x_t)))
sage1Jacobian:{[n;csr;x;lam]
  xt:unitSphere x;
  p:xt*xt;
  ptp:sparsePtMul[n;csr;p];
  sq:sqrt 0f|ptp;
  tgt:unitSphere sq;
  / d(p)/d(x) = 2 * diag(x_t)
  jProb:(2f * xt) */: eye n;
  / d(sqrt(P^T p))/d(p) = diag(1/(2*sqrt)) * P^T
  inv2sq:{$[x>1e-12; 0.5%x; 0f]} each sq;
  ptDense:csrToDenseT[n;csr];
  jSqrt:inv2sq *' ptDense;
  jTgt:jacProj[sq] mmu jSqrt mmu jProb;
  / Geodesic SLERP blend Jacobian on tangent bundle
  ((1f - lam) * eye n) + (lam * jTgt)
 }

/ --- 5. Jacobian Block Recursive Forward Pass ---
/ Recursively propagates (stateBlock X; accumulated Jacobian J_acc; trace table)
/ across hierarchical blocks of (csr; lam; tau)
fwdBlockRec:{[blocks;X;jAcc;depth;trace]
  if[depth >= count blocks; :(X; jAcc; trace)];
  blk:blocks depth;
  n:blk`n; csr:blk`csr; lam:blk`lam; tau:blk`tau;
  / 1) Geometric Riemannian Manifold Attention across block states
  Xattn:riemannAttn[tau; X];
  / 2) SAGE-1 unit-sphere geodesic evolution per state in block
  Xnext:{[n;csr;lam;x]
    xt:unitSphere x;
    tgt:unitSphere sqrt 0f|sparsePtMul[n;csr;xt*xt];
    xn:slerp[xt; tgt; lam];
    if[1e-8 < abs 1f - norm2 xn; '"Unit-sphere invariant violated"];
    xn
  }[n;csr;lam] each Xattn;
  / 3) Block Jacobian composition: J_acc <- J_block * J_acc
  jBlk:sage1Jacobian[n; csr; Xattn 0; lam];
  jNext:jBlk mmu jAcc;
  frob:sqrt sum over jNext*jNext;
  rec:([] depth:enlist depth; mean_norm:enlist avg norm2 each Xnext; jac_frob:enlist frob);
  .z.s[blocks; Xnext; jNext; depth+1; trace,rec]
 }

/ --- 6. Wolfram Flat CSV Pipe Loader (from sage1_geometry_pipe.wl) ---
loadWolframCSV:{[dir]
  csrRaw:("IIF";enlist ",") 0: hsym `$dir,"/sage1_csr.csv";
  pksRaw:("IF";enlist ",") 0: hsym `$dir,"/sage1_peaks.csv";
  x0Raw:("F";enlist ",") 0: hsym `$dir,"/sage1_init.csv";
  csr:update row:row-1, col:col-1 from `row`col`w xcol csrRaw; / Convert 1-based to 0-based (CSV header is row,col,weight)
  n:count x0Raw`x0;
  `n`csr`peaks`x0!(n; csr; pksRaw[`peak_idx]-1; unitSphere x0Raw`x0)
 }

/ --- 7. Standalone Self-Test (5-Vertex Triangulated Patch, 3 Recursive Blocks) ---
runSelfTest:{[]
  n:5;
  csr:([] row:0 0 1 1 1 2 2 3 3 4 4i;
          col:0 1 0 1 2 1 3 2 4 3 4i;
          w:0.5 0.5 0.3 0.4 0.3 0.5 0.5 0.5 0.5 0.5 0.5f);
  x0:unitSphere n#1f;
  x1:unitSphere 0.5 0.3 0.2 0.1 0.1f;
  X0:(x0; x1);
  blocks:(
    `n`csr`lam`tau!(n; csr; 0.30f; 0.50f);
    `n`csr`lam`tau!(n; csr; 0.35f; 0.40f);
    `n`csr`lam`tau!(n; csr; 0.40f; 0.30f)
  );
  res:fwdBlockRec[blocks; X0; "f"$eye n; 0; ()];
  show res 2;
  res
 }

/ Execute self-test when loaded directly
runSelfTest[]
