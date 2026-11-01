Tree-support certificates
=========================

tree_audit.py generates 450 initial chains and 79,113 support records using
exact rational arithmetic. verify_certificates.py reconstructs and replays
every certificate using independent closed-form inverse matrices; it does
not import the generator. Both programs use the Python standard library.

Run from this directory:

  python -B tree_audit.py
  python -B verify_certificates.py

The unified driver at ../../verification/verify_all.py provides these runs
through --regenerate and --replay-tree.

Files
-----
records.jsonl: the 79,113 individually identifiable inputs and certificates.
summary.json: generator summary.
unresolved.json: the unresolved records (empty in the supplied data).
replay_summary.json: independent replay results and SHA-256 of records.jsonl.

The supplied counts are:

  End type    Negative directions    Short path    Failed PD/PSD
  AA                 17115              5877           10083
  BA                  7080              2175            5370
  DA                  7472              2185            5868
  BB                  2400               300            2100
  BD                  2570               310            2400
  DD                  2752               320            2736
  Total              39389             11167           28557

The minimum positive determinant of a two-negative-direction certificate
is 1/32400. These outcomes exclude the prescribed overinclusive local
configurations; they do not assert that each configuration is realizable.

Inputs and reduction assumptions
--------------------------------
The main chain is 0,...,P+ell+Q, with a=P and b=P+ell. Extra endpoint leaves
are not included in these path parameters. End types are AA, BA, DA, BB,
BD, DD. Initial parameters are P,Q=5,...,9 and ell=8,9,10.

The additional marks c,d have no edge between them; an edge would give a
path of length one immediately. The original component C is retained in
full. For c, delete a and generate supports on the resulting finite classical
components; for d, delete b. The function component_supports implements the
prescribed endpoint, near-endpoint, exceptional A-arm, D7/D8 short-leaf,
and all-simple A-cycle support families. It also permits the empty support
and simple endpoint joins across distinct components. Duplicate supports
are removed by their actual edge sets.

The completeness of these support families is a mathematical input. It is
not inferred from agreement with the expected number of records. In the
geometric reduction, diagonal inverse formulas exclude labels at least five;
a double support in one component must be a pure affine A cycle; and Schur
complement conditions restrict joins between distinct components.

When C with a new mark and zero opposite coupling is already indefinite,
zero coupling is allowed as a relaxation. Otherwise the opposite coupling
is required to be nonzero, apart from the prescribed external E7/E8 cases:
a simple leaf near an ordinary original endpoint with near parameter 5 or 6.
Eliminating that external leaf changes the new mark's diagonal in 2G from
2 to 3/2. Both the nonzero-opposite/no-leaf case and the zero-opposite/leaf
case are included. Intersecting or mutually adjacent external leaves, or a
leaf meeting the other mark, supply a short path directly; otherwise the
leaves are eliminated independently.

In the rational congruence, chain and D-leaf diagonals are 2, B-leaf
diagonals are 4, and a heavy-supported new mark has diagonal 4 and support
entry -2. A required opposite coupling is relaxed to magnitude 1. Positive
witness coordinates certify persistence when those opposite couplings are
increased. This is not general monotonicity for arbitrary Gram entries.

Certificate types
-----------------
Eliminate the positive-definite unmarked block S0=C minus {a,b}, retaining
the marks (a,b,c,d) in a 4-by-4 Schur matrix R.

1. Failed necessary definiteness: S0+c and S0+d must be positive definite,
   so R_cc,R_dd>0. The adjacent pairs bc,da,cd must have positive-semidefinite
   2-by-2 principal matrices. A failing pivot or matrix is recorded.
2. Preserved short path: the certificate specifies a path of length at most
   seven in an actual adjacent-pair deletion. The path avoids the other two
   marks and every variable interior segment.
3. Two negative directions: positive rational x,y define the columns
   (1,0,x,0)^T and (0,1,0,y)^T. The restricted 2-by-2 matrix has negative
   first diagonal and positive determinant. Its entries are recorded exactly.

The replayer computes the unmarked inverse independently. For an A trunk
with t nodes, K^{-1}_{ij}=(min(i,j)+1)(t-max(i,j))/(t+1). Eliminating a B/D
decoration gives first diagonal 1 and remaining diagonals 2, with inverse
t-max(i,j). Restoring the B leaf or D leaves gives the remaining entries.

The separate ../vertex_tree_certificate_audit.py also checks the recorded
inequalities and reconstructs every short path without importing the
generator. It does not replace reconstruction of the Schur matrices.

Arbitrary-length coverage
-------------------------
For initial P=9, the variable interior is between nodes 2 and a-2; for
ell=10, between a+2 and b-2; for Q=9, between b+2 and the node two steps
before the right original endpoint. Supports in the stable ranges meet
retained nodes; exceptional supports occur in the smaller fixed cases.

The reduction replaces P,Q,ell by min(P,9),min(Q,9),min(ell,10).
Recovering longer chains inserts nodes only in unmarked simple segments.
The mathematical path-Schur identity makes the retained matrix decrease in
positive-semidefinite order. Further Schur complementation preserves that
order while the eliminated block remains positive definite; loss of that
definiteness already violates the prescribed condition. Thus the failed
definiteness and negative-direction certificates persist. Short-path
certificates explicitly avoid all inserted segments.

This coverage argument is an analytical premise of the finite test. A finite
enumeration alone does not prove the geometric reductions or classification.
