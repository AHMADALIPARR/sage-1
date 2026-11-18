% SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
% Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
%
% Pushdown grammar layer for SAGE-1.
%
% Prolog owns the parser state (a stack plus the previous token) and turns it
% into the `legal` mask for sage1_select/3. Per token the decoder:
%
%   1. steps the sphere state:      sage1_step(Csr, X, Lambda, Y)
%   2. builds the legal mask:        legal_mask(Parser, Mask)
%   3. decodes:                      sage1_select(Y, Mask, selected(Tok, P))
%   4. updates the parser state:     advance(Parser, Tok, Parser1)
%   5. only then conditions:         sage1_condition(Tok, N, X1)
%
% A dead_end from sage1_select makes Prolog backtrack to the previous choice,
% which retries with the next most probable legal token.

:- module(sage1_grammar,
          [ token/2,
            demo_csr/1,
            initial_parser/2,
            legal/2,
            legal_mask/2,
            advance/3,
            decode/6,
            dead_end/1,
            sentence//0
          ]).

:- use_module(sage1).

:- dynamic dead_end/1.

%   States of the demo grammar (shared with rust/sage1 and the Futhark check).
token(0, '<s>').
token(1, '(').
token(2, ')').
token(3, '[').
token(4, ']').
token(5, a).
token(6, '</s>').

states(7).

%   Transition counts of the demo corpus: From-To-Count.
corpus_counts([ 0-3-3, 0-1-2, 0-5-1,
                1-5-4, 1-1-1,
                2-6-3, 2-2-2, 2-4-1,
                3-5-2, 3-1-1,
                4-6-2, 4-5-1,
                5-2-4, 5-6-1,
                6-0-1 ]).

demo_csr(Csr) :-
    states(N),
    corpus_counts(T),
    csr_from_triples(N, T, Csr).

%   Parser state: p(Stack, LastToken, MaxDepth). Stack holds open brackets.
initial_parser(MaxDepth, p([], 0, MaxDepth)).

opener(1, 2).          % ( ... )
opener(3, 4).          % [ ... ]

%!  legal(+Parser, ?Token) is nondet.
legal(p(S, _, Max), Open) :-
    opener(Open, _),
    length(S, D), D < Max.
legal(p([Top|_], _, _), Close) :-
    opener(Top, Close).
legal(p(_, Last, _), 5) :-
    Last =\= 5.                      % no two atoms in a row
legal(p([], Last, _), 6) :-
    Last =\= 0.                      % </s> only on an empty stack, non-empty sentence

legal_mask(P, Mask) :-
    states(N), Last is N - 1,
    numlist(0, Last, Ts),
    maplist(legal_bit(P), Ts, Mask).

legal_bit(P, T, 1) :- legal(P, T), !.
legal_bit(_, _, 0).

%!  advance(+Parser0, +Token, -Parser) is det.  Push/pop the stack.
advance(p(S, _, M), T, p([T|S], T, M)) :- opener(T, _), !.
advance(p([Top|S], _, M), T, p(S, T, M)) :- opener(Top, T), !.
advance(p(S, _, M), T, p(S, T, M)).

%!  decode(+Csr, +Parser0, +Lambda, +MaxLen, -Tokens, -Trace) is semidet.
%
%   Trace is a list of step(Token, Prob, StackAfter, XAfter), one per token,
%   where XAfter is the conditioned sphere state.
decode(Csr, P0, Lambda, MaxLen, Tokens, Trace) :-
    retractall(dead_end(_)),
    csr_states(Csr, N),
    sage1_condition(0, N, X0),
    dfs(Csr, N, Lambda, MaxLen, X0, P0, [], Tokens, Trace).

dfs(Csr, N, Lambda, MaxLen, X, P, Rev, [Tok|Toks], [step(Tok, Prob, S1, X1)|Tr]) :-
    length(Rev, Len), Len < MaxLen,
    sage1_step(Csr, X, Lambda, Y),
    legal_mask(P, Mask),
    choose(Y, Mask, Rev, Tok, Prob),
    advance(P, Tok, P1),                 % parser state first
    P1 = p(S1, _, _),
    sage1_condition(Tok, N, X1),         % then sphere conditioning
    (   Tok =:= 6
    ->  Toks = [], Tr = []
    ;   dfs(Csr, N, Lambda, MaxLen, X1, P1, [Tok|Rev], Toks, Tr)
    ).

choose(Y, Mask, Rev, Tok, Prob) :-
    sage1_select(Y, Mask, R),
    (   R == dead_end
    ->  reverse(Rev, Prefix),
        assertz(dead_end(Prefix)),
        fail
    ;   R = selected(T, P),
        (   Tok = T, Prob = P
        ;   clear_bit(Mask, T, Mask1),
            choose(Y, Mask1, Rev, Tok, Prob)
        )
    ).

clear_bit(Mask, I, Mask1) :-
    nth0(I, Mask, _, Rest),
    nth0(I, Mask1, 0, Rest).

%   Independent DCG for the same language, used to check decoder output.
sentence --> items, [6].
items --> [].
items --> item(Prev), items_after(Prev).
items_after(_) --> [].
items_after(Prev) --> item_ok(Prev, P1), items_after(P1).
item_ok(Prev, P) --> item(P), { \+ (Prev == atom, P == atom) }.
item(atom) --> [5].
item(group) --> [1], inner, [2].
item(group) --> [3], inner, [4].
inner --> [].
inner --> item(P), items_after(P).
