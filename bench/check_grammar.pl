% SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
% Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
%
% Checks the C grammar port and corpus generator used by the perplexity
% benchmark against prolog/sage1_grammar.pl:
%
%   swipl -q -g main -t halt bench/check_grammar.pl -- MASKS CORPUS
%
% MASKS  = output of `sage1_bench --dump-masks D`: every parser state
%          m(Stack, Last, Mask) must equal legal_mask(p(Stack, Last, D), Mask).
% CORPUS = output of `sage1_bench --dump-corpus N Seed D`: every sentence
%          s([0|Tokens]) must parse with sentence//0 and every token must be
%          legal/2 in the parser state reached by advance/3.
% (SAGE1_PL_LIB must point at the built sage1_pl foreign library.)

:- use_module(library(readutil)).
:- prolog_load_context(directory, Dir),
   atom_concat(Dir, '/../prolog', P),
   asserta(user:file_search_path(sage1src, P)).
:- use_module(sage1src(sage1_grammar)).

main :-
    current_prolog_flag(argv, Argv),
    append(_, [Masks, Corpus], Argv),
    check_masks(Masks, NM),
    check_corpus(Corpus, NS, NT),
    format("check_grammar: ~d parser states, ~d sentences (~d tokens) agree with sage1_grammar.pl~n",
           [NM, NS, NT]).

terms(File, Terms) :-
    setup_call_cleanup(open(File, read, In), read_all(In, Terms), close(In)).
read_all(In, Ts) :-
    read_term(In, T, []),
    (   T == end_of_file -> Ts = [] ; Ts = [T|Rest], read_all(In, Rest) ).

check_masks(File, N) :-
    terms(File, [max_depth(D)|Ms]),
    forall(member(m(S, L, Mask), Ms),
           (   legal_mask(p(S, L, D), M1), M1 == Mask -> true
           ;   format(user_error, "mask mismatch at ~w~n", [m(S, L, Mask)]), halt(1) )),
    length(Ms, N), N > 0.

check_corpus(File, NS, NT) :-
    terms(File, [max_depth(D)|Ss]),
    foldl(check_sentence(D), Ss, 0, NT),
    length(Ss, NS), NS > 0.

check_sentence(D, s([0|Toks]), T0, T) :-
    (   phrase(sentence, Toks) -> true
    ;   format(user_error, "not in sentence//0: ~w~n", [Toks]), halt(1) ),
    initial_parser(D, P0),
    (   foldl(legal_step, Toks, P0, _) -> true
    ;   format(user_error, "illegal token in ~w~n", [Toks]), halt(1) ),
    length(Toks, Len),
    T is T0 + Len.

legal_step(Tok, P0, P) :- legal(P0, Tok), advance(P0, Tok, P).
