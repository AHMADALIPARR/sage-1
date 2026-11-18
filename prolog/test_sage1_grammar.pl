% SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
% Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
%
% swipl -g run_tests -t halt prolog/test_sage1_grammar.pl
% (with SAGE1_PL_LIB pointing at the built sage1_pl shared object)

:- use_module(library(plunit)).
:- use_module(sage1).
:- use_module(sage1_grammar).

close_to(A, B) :- abs(A - B) < 1.0e-12.

:- begin_tests(sage1_kernel).

fixture(csr(3, [0,2,3,4], [1,2,2,0], [3,1,1,1])).

test(matches_c_test_sage1) :-
    fixture(M),
    sage1_validate(M, ok),
    sage1_propagate(M, [1,0,0], [_, Q1, Q2]),
    close_to(Q1, 0.75), close_to(Q2, 0.25),
    sage1_step(M, [1,0,0], 1.0, [_, O1, O2]),
    O1sq is O1*O1, O2sq is O2*O2,
    close_to(O1sq, 0.75), close_to(O2sq, 0.25),
    sage1_step(M, [1,0,0], 0.0, [Z0|_]),
    close_to(Z0, 1.0),
    sage1_step(M, [1,0,0], 0.5, H),
    foldl([V,A0,A]>>(A is A0+V*V), H, 0.0, Norm), close_to(Norm, 1.0),
    sage1_select([0.0, O1, O2], [0,1,0], selected(1, P)), close_to(P, 1.0),
    sage1_condition(1, 3, C), C == [0.0, 1.0, 0.0],
    sage1_select(C, [0,0,0], dead_end),
    sage1_validate(csr(3, [0,0,3,4], [1,2,2,0], [3,1,1,1]), bad_matrix).

test(bad_lambda_raises, [throws(error(sage1_error(bad_input), _))]) :-
    fixture(M),
    sage1_step(M, [1,0,0], 2.0, _).

test(geodesic_half) :-
    sage1_geodesic([1,0], [0,1], 0.5, [A,B]),
    S is sqrt(0.5), abs(A-S) < 1.0e-14, abs(B-S) < 1.0e-14.

:- end_tests(sage1_kernel).

:- begin_tests(sage1_grammar).

names(Ts, Ns) :- maplist([T,N]>>token(T,N), Ts, Ns).

test(decode_backtracks_over_dead_end) :-
    demo_csr(Csr),
    initial_parser(4, P0),
    once(decode(Csr, P0, 1.0, 16, Toks, Trace)),
    names(Toks, Names),
    Names == ['[', '(', a, ')', ']', '</s>'],
    findall(D, dead_end(D), Dead),
    Dead == [[3, 5]],
    maplist([step(_, Pr, _, _), Pr]>>true, Trace, Probs),
    maplist(close_to, Probs, [0.5, 1.0, 0.8, 1.0, 1.0, 0.6666666666666666]).

test(condition_follows_parser_update) :-
    demo_csr(Csr),
    initial_parser(4, P0),
    once(decode(Csr, P0, 1.0, 16, _, Trace)),
    maplist([step(_, _, S, _), S]>>true, Trace, Stacks),
    Stacks == [[3], [1,3], [1,3], [3], [], []],
    forall(member(step(T, _, _, X), Trace),
           ( nth0(T, X, 1.0), sum_list(X, 1.0) )).

test(select_rejects_dead_end) :-
    demo_csr(Csr),
    initial_parser(4, P0),
    advance(P0, 3, P1), advance(P1, 5, P2),
    sage1_condition(5, 7, X),
    sage1_step(Csr, X, 1.0, Y),
    legal_mask(P2, Mask),
    Mask == [0,1,0,1,1,0,0],
    sage1_select(Y, Mask, dead_end).

test(decoded_sentences_are_in_the_language) :-
    demo_csr(Csr),
    forall(( between(1, 6, Depth), member(L, [1.0, 0.9, 0.75]),
             initial_parser(Depth, P0),
             once(decode(Csr, P0, L, 32, Toks, _)) ),
           phrase(sentence, Toks)).

test(dcg_rejects_unbalanced, [fail]) :-
    phrase(sentence, [3, 1, 5, 4, 2, 6]).

:- end_tests(sage1_grammar).
