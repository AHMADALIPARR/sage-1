-- SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
-- Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
--
-- SAGE-1 array backend in Futhark: data-parallel versions of the C11 kernel
-- in native/sage1. Every entry returns a status code first, using the same
-- values as sage1_status (0 ok, 1 bad input, 2 bad matrix, 3 dead end,
-- 4 numeric). CSR column indices are i64 here, and the matrix is expected to
-- have passed sage1_validate (rows nonempty, positive counts).
--
--   propagate : q = P^T p                       (segmented CSR scatter-reduce)
--   step      : SLERP(x, sqrt(P^T (x*x)), λ)
--   geodesic  : SLERP(x, y, λ) on the nonnegative unit sphere
--   select    : legal argmax of x*x, ties -> smallest index
--   condition : one-hot collapse onto the selected state

def ok: i32 = 0
def bad_input: i32 = 1
def bad_matrix: i32 = 2
def dead_end: i32 = 3
def numeric: i32 = 4

def isfinite (a: f64): bool = !(f64.isnan a || f64.isinf a)

def finite_unit [n] (v: [n]f64): bool =
  n > 0
  && all (\a -> isfinite a && a >= 0) v
  && (let s = f64.sum (map (\a -> a * a) v)
      in isfinite s && f64.abs (s - 1) <= 1e-8)

-- Source row of every edge, from the CSR row pointer.
def edge_rows [n1] (e: i64) (row_ptr: [n1]i64): [e]i64 =
  let n = n1 - 1
  let starts = map (\i -> if i < e then i else -1) row_ptr[1:n]
  let flags = scatter (replicate e 0i64) starts (replicate (n - 1) 1i64)
  in scan (+) 0 flags

def csr_shape_ok [n1][e] (n: i64) (row_ptr: [n1]i64) (col: [e]i64) (counts: [e]u64): bool =
  n1 == n + 1 && e > 0 && row_ptr[0] == 0 && row_ptr[n] == e
  && all (\i -> row_ptr[i] < row_ptr[i + 1]) (iota n)
  && all (\c -> c >= 0 && c < n) col
  && all (> 0) counts

def propagate_impl [n1][e][n] (row_ptr: [n1]i64) (col: [e]i64) (counts: [e]u64) (p: [n]f64): (i32, [n]f64) =
  if !(csr_shape_ok n row_ptr col counts) then (bad_matrix, replicate n 0)
  else if !(all (\a -> isfinite a && a >= 0) p) then (bad_input, replicate n 0)
  else
    let mass = f64.sum p
    in if !(isfinite mass) || f64.abs (mass - 1) > 1e-8 then (bad_input, replicate n 0)
    else
      let rows = edge_rows e row_ptr
      let totals = hist (+) 0u64 n rows counts
      let flow = map2 (\r k -> p[r] * (f64.u64 k / f64.u64 totals[r])) rows counts
      let q = hist (+) 0f64 n col flow
      let sum = f64.sum q
      in if !(sum > 0) || !(isfinite sum) then (numeric, replicate n 0)
      else (ok, map (/ sum) q)

def geodesic_impl [n] (x: [n]f64) (y: [n]f64) (lambda: f64): (i32, [n]f64) =
  if !(finite_unit x) || !(finite_unit y) || !(isfinite lambda) || lambda < 0 || lambda > 1
  then (bad_input, replicate n 0)
  else
    let dot = f64.max (-1) (f64.min 1 (f64.sum (map2 (*) x y)))
    let theta = f64.acos dot
    let sin_theta = f64.sin theta
    let (a, b) = if theta < 1e-7 then (1 - lambda, lambda)
                 else (f64.sin ((1 - lambda) * theta) / sin_theta,
                       f64.sin (lambda * theta) / sin_theta)
    let out = map2 (\xi yi -> let o = a * xi + b * yi
                              in if o < 0 && o > -1e-14 then 0 else o) x y
    let norm2 = f64.sum (map (\o -> o * o) out)
    in if any (\o -> !(isfinite o) || o < 0) out || !(norm2 > 0) || !(isfinite norm2)
       then (numeric, replicate n 0)
       else let inv = 1 / f64.sqrt norm2 in (ok, map (* inv) out)

entry propagate [n1][e][n] (row_ptr: [n1]i64) (col: [e]i64) (counts: [e]u64) (p: [n]f64): (i32, [n]f64) =
  propagate_impl row_ptr col counts p

entry geodesic [n] (x: [n]f64) (y: [n]f64) (lambda: f64): (i32, [n]f64) =
  geodesic_impl x y lambda

entry step [n1][e][n] (row_ptr: [n1]i64) (col: [e]i64) (counts: [e]u64) (x: [n]f64) (lambda: f64): (i32, [n]f64) =
  if !(finite_unit x) then (bad_input, replicate n 0)
  else
    let (s, q) = propagate_impl row_ptr col counts (map (\a -> a * a) x)
    in if s != ok then (s, replicate n 0)
       else geodesic_impl x (map f64.sqrt q) lambda

-- Fused multi-step evolution: x_{t+1} = step x_t, for `steps` iterations.
entry evolve [n1][e][n] (row_ptr: [n1]i64) (col: [e]i64) (counts: [e]u64) (x0: [n]f64) (lambda: f64) (steps: i64): (i32, [n]f64) =
  loop (s, x) = (ok, x0) for _i < steps do
    if s != ok then (s, x)
    else if !(finite_unit x) then (bad_input, x)
    else
      let (s1, q) = propagate_impl row_ptr col counts (map (\a -> a * a) x)
      in if s1 != ok then (s1, x)
         else geodesic_impl x (map f64.sqrt q) lambda

entry select [n] (x: [n]f64) (legal: [n]u8): (i32, i64, f64) =
  if !(finite_unit x) then (bad_input, 0, 0)
  else
    let w = map2 (\xi l -> if l != 0 then xi * xi else -1) x legal
    let total = f64.sum (map (\a -> if a >= 0 then a else 0) w)
    let pick (wa: f64, ia: i64) (wb: f64, ib: i64) =
      if wa > wb then (wa, ia) else if wb > wa then (wb, ib)
      else (wa, i64.min ia ib)
    let (best, idx) = reduce_comm pick (-f64.inf, i64.highest) (zip w (iota n))
    in if !(total > 0) then (dead_end, 0, 0)
       else (ok, idx, best / total)

entry condition (selected: i64) (n: i64): (i32, []f64) =
  if n <= 0 || selected < 0 || selected >= n then (bad_input, replicate (i64.max n 0) 0)
  else (ok, map (\i -> if i == selected then 1 else 0) (iota n))
