// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
//
// Mirrors every assertion of native/sage1/test_sage1.c through the safe API.

use sage1::{condition, geodesic, select, Csr, Error};

fn fixture() -> Csr {
    Csr::new(vec![0, 2, 3, 4], vec![1, 2, 2, 0], vec![3, 1, 1, 1]).expect("valid CSR")
}

#[test]
fn matches_c_test_sage1() {
    let m = fixture();
    assert_eq!(m.validate(), Ok(()));

    let q = m.propagate(&[1.0, 0.0, 0.0]).unwrap();
    assert!((q[1] - 0.75).abs() < 1e-12 && (q[2] - 0.25).abs() < 1e-12);

    let x = [1.0, 0.0, 0.0];
    let out = m.step(&x, 1.0).unwrap();
    assert!((out[1] * out[1] - 0.75).abs() < 1e-12);
    assert!((out[2] * out[2] - 0.25).abs() < 1e-12);
    let target = out.clone();

    let out = m.step(&x, 0.0).unwrap();
    assert!((out[0] - 1.0).abs() < 1e-12);

    let out = m.step(&x, 0.5).unwrap();
    assert!((out.iter().map(|v| v * v).sum::<f64>() - 1.0).abs() < 1e-12);

    let sel = select(&target, &[false, true, false]).unwrap();
    assert_eq!(sel.index, 1);
    assert!((sel.probability - 1.0).abs() < 1e-12);

    let out = condition(sel.index, 3).unwrap();
    assert_eq!(out, vec![0.0, 1.0, 0.0]);

    assert_eq!(select(&out, &[false, false, false]), Err(Error::DeadEnd));

    let invalid = Csr::new(vec![0, 0, 3, 4], vec![1, 2, 2, 0], vec![3, 1, 1, 1]);
    assert_eq!(invalid, Err(Error::BadMatrix));
}

#[test]
fn rejects_bad_shapes_before_ffi() {
    let m = fixture();
    assert!(matches!(m.propagate(&[1.0, 0.0]), Err(Error::LengthMismatch { expected: 3, found: 2 })));
    assert!(matches!(m.step(&[1.0, 0.0, 0.0, 0.0], 0.5), Err(Error::LengthMismatch { .. })));
    assert!(matches!(select(&[1.0], &[true, false]), Err(Error::LengthMismatch { .. })));
    assert_eq!(Csr::new(vec![0, 2, 3, 9], vec![1, 2, 2, 0], vec![3, 1, 1, 1]), Err(Error::BadInput));
    assert_eq!(Csr::new(vec![0, 7, 3, 4], vec![1, 2, 2, 0], vec![3, 1, 1, 1]), Err(Error::BadMatrix));
    assert_eq!(Csr::new(vec![0, 2, 3, 4], vec![1, 3, 2, 0], vec![3, 1, 1, 1]), Err(Error::BadMatrix));
    assert_eq!(Csr::new(vec![0, 2, 3, 4], vec![1, 2, 2, 0], vec![3, 0, 1, 1]), Err(Error::BadMatrix));
}

#[test]
fn kernel_errors_map_to_rust() {
    let m = fixture();
    assert_eq!(m.step(&[1.0, 0.0, 0.0], 1.5), Err(Error::BadInput));
    assert_eq!(m.step(&[0.5, 0.5, 0.0], 0.5), Err(Error::BadInput));
    assert_eq!(m.propagate(&[0.5, 0.0, 0.0]), Err(Error::BadInput));
    assert_eq!(condition(3, 3), Err(Error::BadInput));
}

#[test]
fn geodesic_quarter_arc() {
    let a = [1.0, 0.0];
    let b = [0.0, 1.0];
    let g = geodesic(&a, &b, 0.25).unwrap();
    let t = std::f64::consts::PI / 8.0;
    assert!((g[0] - t.cos()).abs() < 1e-14 && (g[1] - t.sin()).abs() < 1e-14);
}

#[test]
fn lambda_one_is_markov_over_many_steps() {
    let m = fixture();
    let mut x = vec![1.0, 0.0, 0.0];
    let mut p = vec![1.0, 0.0, 0.0];
    for _ in 0..50 {
        x = m.step(&x, 1.0).unwrap();
        p = m.propagate(&p).unwrap();
        for i in 0..3 {
            assert!((x[i] * x[i] - p[i]).abs() < 1e-12);
        }
    }
}

#[test]
fn from_triples_builds_same_csr() {
    let t = Csr::from_triples(3, &[(0, 1, 3), (0, 2, 1), (1, 2, 1), (2, 0, 1)]).unwrap();
    assert_eq!(t, fixture());
}
