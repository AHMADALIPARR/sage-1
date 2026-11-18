% SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
% Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
%
% SWI-Prolog binding for the SAGE-1 C11 kernel (foreign library sage1_pl).
% Set SAGE1_PL_LIB to the built shared object, or put it on the
% foreign/1 search path.

:- module(sage1,
          [ sage1_validate/2,
            sage1_propagate/3,
            sage1_step/4,
            sage1_geodesic/4,
            sage1_select/3,
            sage1_condition/3,
            csr_from_triples/3,
            csr_states/2
          ]).

load_sage1_foreign :-
    (   getenv('SAGE1_PL_LIB', Lib)
    ->  use_foreign_library(Lib)
    ;   use_foreign_library(foreign(sage1_pl))
    ).

:- initialization(load_sage1_foreign, now).

%!  csr_from_triples(+States, +Triples, -Csr) is det.
%
%   Build csr(States, RowPtr, ColIdx, Counts) from From-To-Count triples,
%   grouping edges by source row in order of appearance.
csr_from_triples(N, Triples, csr(N, RowPtr, Cols, Counts)) :-
    Last is N - 1,
    numlist(0, Last, Rows),
    foldl(row_edges(Triples), Rows, Edges0, []),
    Edges0 = Edges,
    pairs_cols_counts(Edges, Cols, Counts),
    foldl(row_ptr(Triples), Rows, Ptrs, 0, _),
    append([0], Ptrs, RowPtr).

row_edges(Triples, Row, Edges0, Edges) :-
    findall(To-C, member(Row-To-C, Triples), RowEdges),
    append(RowEdges, Edges, Edges0).

row_ptr(Triples, Row, Ptr, Acc0, Ptr) :-
    aggregate_all(count, member(Row-_-_, Triples), K),
    Ptr is Acc0 + K.

pairs_cols_counts([], [], []).
pairs_cols_counts([T-C|Es], [T|Ts], [C|Cs]) :- pairs_cols_counts(Es, Ts, Cs).

csr_states(csr(N, _, _, _), N).
