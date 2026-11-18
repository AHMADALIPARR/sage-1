// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab

//! Safe Rust runtime for the SAGE-1 C11 kernel.
//!
//! SAGE-1 evolves a state `x` on the unit sphere with
//! `x_{t+1} = SLERP(x_t, sqrt(P^T (x_t ⊙ x_t)), λ)`, where `P` is the
//! row-normalized sparse count matrix of grammatical transitions. A grammar
//! supplies legal-token masks, [`select`] decodes deterministically, and
//! [`condition`] collapses the state onto the chosen branch after the parser
//! state has been updated. [`decode`] runs that loop with backtracking over
//! dead ends for any [`Grammar`].
//!
//! Every slice length is checked in Rust before crossing the FFI boundary, and
//! the scratch/output buffers handed to C are always distinct allocations, so
//! the C aliasing rules can never be violated from safe code.

pub mod ffi;
pub mod grammar;

use std::fmt;

/// Errors reported by the kernel (`sage1_status`) or by the Rust wrapper.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Error {
    /// `SAGE1_BAD_INPUT`: non-unit, negative or non-finite vectors, bad λ, bad sizes.
    BadInput,
    /// `SAGE1_BAD_MATRIX`: empty rows, column out of range, zero counts, overflow.
    BadMatrix,
    /// `SAGE1_DEAD_END`: no probability mass on any legal transition.
    DeadEnd,
    /// `SAGE1_NUMERIC`: the propagation or geodesic left the valid domain.
    Numeric,
    /// A slice length did not match the number of states (checked in Rust).
    LengthMismatch { expected: usize, found: usize },
    /// An unknown status code came back from C.
    Unknown(i32),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::BadInput => write!(f, "SAGE-1: bad input"),
            Error::BadMatrix => write!(f, "SAGE-1: bad CSR matrix"),
            Error::DeadEnd => write!(f, "SAGE-1: dead end (no legal transition has mass)"),
            Error::Numeric => write!(f, "SAGE-1: numeric failure"),
            Error::LengthMismatch { expected, found } => {
                write!(f, "SAGE-1: expected length {}, found {}", expected, found)
            }
            Error::Unknown(c) => write!(f, "SAGE-1: unknown status {}", c),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

fn check(status: ffi::sage1_status) -> Result<()> {
    match status {
        ffi::SAGE1_OK => Ok(()),
        ffi::SAGE1_BAD_INPUT => Err(Error::BadInput),
        ffi::SAGE1_BAD_MATRIX => Err(Error::BadMatrix),
        ffi::SAGE1_DEAD_END => Err(Error::DeadEnd),
        ffi::SAGE1_NUMERIC => Err(Error::Numeric),
        other => Err(Error::Unknown(other)),
    }
}

fn expect_len(expected: usize, found: usize) -> Result<()> {
    if expected == found {
        Ok(())
    } else {
        Err(Error::LengthMismatch { expected, found })
    }
}

/// An owned, validated CSR transition-count matrix. Row `i` holds the
/// transitions `i -> col_idx[k]` with positive integer counts.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Csr {
    row_ptr: Vec<usize>,
    col_idx: Vec<u32>,
    counts: Vec<u64>,
}

impl Csr {
    /// Build from raw CSR arrays; validated by `sage1_validate`.
    pub fn new(row_ptr: Vec<usize>, col_idx: Vec<u32>, counts: Vec<u64>) -> Result<Self> {
        if row_ptr.len() < 2 || col_idx.len() != counts.len() {
            return Err(Error::BadInput);
        }
        if row_ptr.last() != Some(&col_idx.len()) {
            return Err(Error::BadInput);
        }
        let m = Csr { row_ptr, col_idx, counts };
        // row_ptr entries may still point past the edge arrays; the C validator
        // checks a<b<=edges per row, so make sure no index escapes first.
        if m.row_ptr.iter().any(|&r| r > m.col_idx.len()) {
            return Err(Error::BadMatrix);
        }
        m.validate()?;
        Ok(m)
    }

    /// Build from `(from, to, count)` triples over `states` states.
    /// Triples are grouped by source row in order of appearance.
    pub fn from_triples(states: usize, triples: &[(u32, u32, u64)]) -> Result<Self> {
        let mut rows: Vec<Vec<(u32, u64)>> = vec![Vec::new(); states];
        for &(from, to, count) in triples {
            let r = rows.get_mut(from as usize).ok_or(Error::BadMatrix)?;
            r.push((to, count));
        }
        let mut row_ptr = Vec::with_capacity(states + 1);
        let mut col_idx = Vec::with_capacity(triples.len());
        let mut counts = Vec::with_capacity(triples.len());
        row_ptr.push(0);
        for r in rows {
            for (to, c) in r {
                col_idx.push(to);
                counts.push(c);
            }
            row_ptr.push(col_idx.len());
        }
        Csr::new(row_ptr, col_idx, counts)
    }

    pub fn states(&self) -> usize {
        self.row_ptr.len() - 1
    }
    pub fn edges(&self) -> usize {
        self.col_idx.len()
    }
    pub fn row_ptr(&self) -> &[usize] {
        &self.row_ptr
    }
    pub fn col_idx(&self) -> &[u32] {
        &self.col_idx
    }
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }

    fn raw(&self) -> ffi::sage1_csr {
        ffi::sage1_csr {
            states: self.states(),
            edges: self.edges(),
            row_ptr: self.row_ptr.as_ptr(),
            col_idx: self.col_idx.as_ptr(),
            counts: self.counts.as_ptr(),
        }
    }

    /// `sage1_validate`.
    pub fn validate(&self) -> Result<()> {
        let raw = self.raw();
        check(unsafe { ffi::sage1_validate(&raw) })
    }

    /// `q = P^T p` for a probability vector `p` (`sage1_propagate`).
    pub fn propagate(&self, p: &[f64]) -> Result<Vec<f64>> {
        let n = self.states();
        expect_len(n, p.len())?;
        let mut q = vec![0.0; n];
        let raw = self.raw();
        check(unsafe { ffi::sage1_propagate(&raw, p.as_ptr(), q.as_mut_ptr()) })?;
        Ok(q)
    }

    /// One SAGE-1 step `SLERP(x, sqrt(P^T (x⊙x)), λ)` (`sage1_step`).
    pub fn step(&self, x: &[f64], lambda: f64) -> Result<Vec<f64>> {
        let mut out = vec![0.0; self.states()];
        let mut scratch = Scratch::new(self.states());
        self.step_into(x, lambda, &mut out, &mut scratch)?;
        Ok(out)
    }

    /// Allocation-free [`Csr::step`] using caller-owned output and scratch.
    pub fn step_into(&self, x: &[f64], lambda: f64, out: &mut [f64], scratch: &mut Scratch) -> Result<()> {
        let n = self.states();
        expect_len(n, x.len())?;
        expect_len(n, out.len())?;
        expect_len(n, scratch.q.len())?;
        expect_len(n, scratch.y.len())?;
        let raw = self.raw();
        // `x` is a shared borrow and `out`/scratch are unique borrows of
        // distinct buffers, so none of them alias.
        check(unsafe {
            ffi::sage1_step(
                &raw,
                x.as_ptr(),
                lambda,
                out.as_mut_ptr(),
                scratch.q.as_mut_ptr(),
                scratch.y.as_mut_ptr(),
            )
        })
    }
}

/// Scratch buffers for [`Csr::step_into`].
#[derive(Debug, Clone)]
pub struct Scratch {
    q: Vec<f64>,
    y: Vec<f64>,
}

impl Scratch {
    pub fn new(states: usize) -> Self {
        Scratch { q: vec![0.0; states], y: vec![0.0; states] }
    }
}

/// Spherical interpolation between unit vectors in the nonnegative orthant
/// (`sage1_geodesic`).
pub fn geodesic(x: &[f64], y: &[f64], lambda: f64) -> Result<Vec<f64>> {
    expect_len(x.len(), y.len())?;
    let mut out = vec![0.0; x.len()];
    check(unsafe { ffi::sage1_geodesic(x.as_ptr(), y.as_ptr(), x.len(), lambda, out.as_mut_ptr()) })?;
    Ok(out)
}

/// Result of decoding a state against a legal mask.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Selection {
    /// Selected state (argmax of `x_i^2` over legal states, ties → smallest index).
    pub index: usize,
    /// Its probability renormalized over the legal states.
    pub probability: f64,
}

/// Deterministic spherical decoding restricted to `legal` (`sage1_select`).
/// Returns [`Error::DeadEnd`] when no legal state carries probability mass.
pub fn select(x: &[f64], legal: &[bool]) -> Result<Selection> {
    expect_len(x.len(), legal.len())?;
    let mask: Vec<u8> = legal.iter().map(|&b| b as u8).collect();
    let mut index = 0usize;
    let mut probability = 0.0f64;
    check(unsafe { ffi::sage1_select(x.as_ptr(), mask.as_ptr(), x.len(), &mut index, &mut probability) })?;
    Ok(Selection { index, probability })
}

/// One-hot state collapsed onto `selected` (`sage1_condition`).
pub fn condition(selected: usize, n: usize) -> Result<Vec<f64>> {
    let mut x = vec![0.0; n];
    condition_in_place(selected, &mut x)?;
    Ok(x)
}

/// In-place [`condition`].
pub fn condition_in_place(selected: usize, x: &mut [f64]) -> Result<()> {
    check(unsafe { ffi::sage1_condition(selected, x.len(), x.as_mut_ptr()) })
}

pub use grammar::{decode, Decoded, Grammar};
