// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
//
// Raw declarations matching native/sage1/sage1.h.

#![allow(non_camel_case_types)]

use std::os::raw::{c_double, c_int};

#[repr(C)]
pub struct sage1_csr {
    pub states: usize,
    pub edges: usize,
    pub row_ptr: *const usize,
    pub col_idx: *const u32,
    pub counts: *const u64,
}

/// `sage1_status` is a plain C enum, i.e. an `int` on every supported ABI.
pub type sage1_status = c_int;
pub const SAGE1_OK: sage1_status = 0;
pub const SAGE1_BAD_INPUT: sage1_status = 1;
pub const SAGE1_BAD_MATRIX: sage1_status = 2;
pub const SAGE1_DEAD_END: sage1_status = 3;
pub const SAGE1_NUMERIC: sage1_status = 4;

extern "C" {
    pub fn sage1_validate(m: *const sage1_csr) -> sage1_status;
    pub fn sage1_propagate(m: *const sage1_csr, p: *const c_double, q: *mut c_double) -> sage1_status;
    pub fn sage1_geodesic(
        x: *const c_double,
        y: *const c_double,
        n: usize,
        lambda: c_double,
        out: *mut c_double,
    ) -> sage1_status;
    pub fn sage1_step(
        m: *const sage1_csr,
        x: *const c_double,
        lambda: c_double,
        out: *mut c_double,
        q: *mut c_double,
        y: *mut c_double,
    ) -> sage1_status;
    pub fn sage1_select(
        x: *const c_double,
        legal: *const u8,
        n: usize,
        selected: *mut usize,
        probability: *mut c_double,
    ) -> sage1_status;
    pub fn sage1_condition(selected: usize, n: usize, x: *mut c_double) -> sage1_status;
}
