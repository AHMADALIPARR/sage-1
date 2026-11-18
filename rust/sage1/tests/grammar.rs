// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
//
// Grammar-constrained decoding with the demo pushdown grammar. The expected
// trace is the same one the Prolog layer produces (prolog/test_sage1_grammar.pl).

use sage1::grammar::brackets::{self, Brackets};
use sage1::{decode, select, Error, Grammar};

#[test]
fn decodes_balanced_sentence_with_backtracking() {
    let csr = brackets::csr().unwrap();
    let d = decode(&csr, &Brackets::new(4), brackets::BOS, 1.0, 16).unwrap();
    let words: Vec<&str> = d.tokens.iter().map(|&t| brackets::NAMES[t]).collect();
    assert_eq!(words, ["[", "(", "a", ")", "]", "</s>"]);
    // "[ a" is rejected: after an atom inside "[", the corpus only continues
    // with ")" or "</s>", both illegal, so sage1_select reports a dead end.
    assert_eq!(d.dead_ends, vec![vec![brackets::LBRACK, brackets::ATOM]]);
    let expect = [0.5, 1.0, 0.8, 1.0, 1.0, 2.0 / 3.0];
    for (p, e) in d.probabilities.iter().zip(expect) {
        assert!((p - e).abs() < 1e-12, "{} vs {}", p, e);
    }
}

#[test]
fn select_rejects_the_dead_end_directly() {
    let csr = brackets::csr().unwrap();
    let mut g = Brackets::new(4);
    g.advance(brackets::LBRACK);
    g.advance(brackets::ATOM);
    let x = sage1::condition(brackets::ATOM, brackets::STATES).unwrap();
    let y = csr.step(&x, 1.0).unwrap();
    let mut mask = vec![false; brackets::STATES];
    g.legal(&mut mask);
    assert!(mask[brackets::RBRACK] && !mask[brackets::RPAREN] && !mask[brackets::EOS]);
    assert_eq!(select(&y, &mask), Err(Error::DeadEnd));
}

#[test]
fn every_decoded_sentence_is_balanced() {
    let csr = brackets::csr().unwrap();
    for depth in 1..=6 {
        for &lambda in &[1.0, 0.9, 0.75] {
            if let Ok(d) = decode(&csr, &Brackets::new(depth), brackets::BOS, lambda, 32) {
                let mut g = Brackets::new(depth);
                for &t in &d.tokens {
                    let mut mask = vec![false; brackets::STATES];
                    g.legal(&mut mask);
                    assert!(mask[t]);
                    g.advance(t);
                }
                assert!(g.stack.is_empty());
                assert_eq!(*d.tokens.last().unwrap(), brackets::EOS);
            }
        }
    }
}
