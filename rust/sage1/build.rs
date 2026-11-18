// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
//
// Compiles native/sage1/sage1.c into a static library with the system C
// compiler (no build dependencies). Honors $CC, $AR and $SAGE1_CFLAGS.

use std::env;
use std::path::PathBuf;
use std::process::Command;

fn run(cmd: &mut Command) {
    let status = cmd
        .status()
        .unwrap_or_else(|e| panic!("failed to spawn {:?}: {}", cmd, e));
    assert!(status.success(), "command failed: {:?}", cmd);
}

fn main() {
    let manifest = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let native = manifest.join("../../native/sage1").canonicalize().expect("native/sage1 not found");
    let src = native.join("sage1.c");
    let out = PathBuf::from(env::var("OUT_DIR").unwrap());
    let obj = out.join("sage1.o");
    let lib = out.join("libsage1.a");

    let cc = env::var("CC").unwrap_or_else(|_| "cc".into());
    let ar = env::var("AR").unwrap_or_else(|_| "ar".into());

    let mut compile = Command::new(&cc);
    compile.args(["-std=c11", "-O2", "-fPIC", "-Wall", "-Wextra", "-Werror", "-pedantic", "-c"]);
    if let Ok(extra) = env::var("SAGE1_CFLAGS") {
        compile.args(extra.split_whitespace());
    }
    compile.arg(&src).arg("-o").arg(&obj);
    run(&mut compile);

    let _ = std::fs::remove_file(&lib);
    run(Command::new(&ar).arg("crs").arg(&lib).arg(&obj));

    println!("cargo:rustc-link-search=native={}", out.display());
    println!("cargo:rustc-link-lib=static=sage1");
    println!("cargo:rustc-link-lib=m");
    println!("cargo:include={}", native.display());
    println!("cargo:rerun-if-changed={}", src.display());
    println!("cargo:rerun-if-changed={}", native.join("sage1.h").display());
    println!("cargo:rerun-if-env-changed=CC");
    println!("cargo:rerun-if-env-changed=AR");
    println!("cargo:rerun-if-env-changed=SAGE1_CFLAGS");
}
