// SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
// Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab

//! Grammar-constrained decoding on top of the kernel.
//!
//! The loop per token is: `step` the sphere state, ask the grammar for the
//! legal mask, `select`, update the parser state with the chosen token, and
//! only then `condition` the sphere state onto it. A `DeadEnd` from `select`
//! makes the decoder backtrack to the previous choice and try the next legal
//! token, in decreasing probability order.

use crate::{condition, select, Csr, Error, Result, Scratch};

/// A parser whose state decides which next states are legal.
pub trait Grammar: Clone {
    /// Fill `mask[i] = true` for every state that may come next.
    fn legal(&self, mask: &mut [bool]);
    /// Update the parser state (e.g. push/pop the stack) after `token`.
    fn advance(&mut self, token: usize);
    /// True if decoding stops after `token` (an end-of-sequence state).
    fn is_final(&self, token: usize) -> bool;
}

/// Output of [`decode`].
#[derive(Debug, Clone, PartialEq)]
pub struct Decoded {
    /// Emitted states, ending with a final state.
    pub tokens: Vec<usize>,
    /// Probability of each emitted state among the legal ones at its step.
    pub probabilities: Vec<f64>,
    /// Prefixes at which `select` reported a dead end.
    pub dead_ends: Vec<Vec<usize>>,
}

struct Search<'a, G: Grammar> {
    csr: &'a Csr,
    lambda: f64,
    max_len: usize,
    scratch: Scratch,
    tokens: Vec<usize>,
    probs: Vec<f64>,
    dead_ends: Vec<Vec<usize>>,
    _g: std::marker::PhantomData<G>,
}

impl<'a, G: Grammar> Search<'a, G> {
    fn dfs(&mut self, x: &[f64], g: &G) -> Result<bool> {
        if self.tokens.len() >= self.max_len {
            return Ok(false);
        }
        let n = self.csr.states();
        let mut y = vec![0.0; n];
        self.csr.step_into(x, self.lambda, &mut y, &mut self.scratch)?;
        let mut mask = vec![false; n];
        g.legal(&mut mask);
        loop {
            let sel = match select(&y, &mask) {
                Ok(s) => s,
                Err(Error::DeadEnd) => {
                    self.dead_ends.push(self.tokens.clone());
                    return Ok(false);
                }
                Err(e) => return Err(e),
            };
            let tok = sel.index;
            // Parser state first, then sphere conditioning.
            let mut g2 = g.clone();
            g2.advance(tok);
            let x2 = condition(tok, n)?;
            self.tokens.push(tok);
            self.probs.push(sel.probability);
            if g2.is_final(tok) || self.dfs(&x2, &g2)? {
                return Ok(true);
            }
            self.tokens.pop();
            self.probs.pop();
            mask[tok] = false;
        }
    }
}

/// Decode from `start` (a one-hot state) under `grammar`, at most `max_len`
/// tokens. Returns `Err(DeadEnd)` if every branch dead-ends.
pub fn decode<G: Grammar>(csr: &Csr, grammar: &G, start: usize, lambda: f64, max_len: usize) -> Result<Decoded> {
    let x0 = condition(start, csr.states())?;
    let mut s: Search<G> = Search {
        csr,
        lambda,
        max_len,
        scratch: Scratch::new(csr.states()),
        tokens: Vec::new(),
        probs: Vec::new(),
        dead_ends: Vec::new(),
        _g: std::marker::PhantomData,
    };
    if s.dfs(&x0, grammar)? {
        Ok(Decoded { tokens: s.tokens, probabilities: s.probs, dead_ends: s.dead_ends })
    } else {
        Err(Error::DeadEnd)
    }
}

/// The demo pushdown grammar shared with the Prolog layer
/// (`prolog/sage1_grammar.pl`): balanced `( )` / `[ ]` with atoms `a`,
/// no two atoms in a row, bounded depth, and `</s>` only on an empty stack.
pub mod brackets {
    use super::Grammar;
    use crate::{Csr, Result};

    pub const BOS: usize = 0;
    pub const LPAREN: usize = 1;
    pub const RPAREN: usize = 2;
    pub const LBRACK: usize = 3;
    pub const RBRACK: usize = 4;
    pub const ATOM: usize = 5;
    pub const EOS: usize = 6;
    pub const STATES: usize = 7;
    pub const NAMES: [&str; STATES] = ["<s>", "(", ")", "[", "]", "a", "</s>"];

    /// Transition counts of the demo corpus (identical in Prolog and Futhark tests).
    pub const COUNTS: [(u32, u32, u64); 15] = [
        (0, 3, 3), (0, 1, 2), (0, 5, 1),
        (1, 5, 4), (1, 1, 1),
        (2, 6, 3), (2, 2, 2), (2, 4, 1),
        (3, 5, 2), (3, 1, 1),
        (4, 6, 2), (4, 5, 1),
        (5, 2, 4), (5, 6, 1),
        (6, 0, 1),
    ];

    pub fn csr() -> Result<Csr> {
        Csr::from_triples(STATES, &COUNTS)
    }

    #[derive(Debug, Clone, PartialEq, Eq)]
    pub struct Brackets {
        pub stack: Vec<usize>,
        pub last: usize,
        pub max_depth: usize,
    }

    impl Brackets {
        pub fn new(max_depth: usize) -> Self {
            Brackets { stack: Vec::new(), last: BOS, max_depth }
        }
    }

    impl Grammar for Brackets {
        fn legal(&self, mask: &mut [bool]) {
            for m in mask.iter_mut() {
                *m = false;
            }
            let open_ok = self.stack.len() < self.max_depth;
            mask[LPAREN] = open_ok;
            mask[LBRACK] = open_ok;
            mask[ATOM] = self.last != ATOM;
            match self.stack.last() {
                Some(&LPAREN) => mask[RPAREN] = true,
                Some(&LBRACK) => mask[RBRACK] = true,
                _ => {}
            }
            mask[EOS] = self.stack.is_empty() && self.last != BOS;
        }
        fn advance(&mut self, token: usize) {
            match token {
                LPAREN | LBRACK => self.stack.push(token),
                RPAREN | RBRACK => {
                    self.stack.pop();
                }
                _ => {}
            }
            self.last = token;
        }
        fn is_final(&self, token: usize) -> bool {
            token == EOS
        }
    }
}
