# Phase E: the Sangaku CAS bridge — verified against the exact sources

Updated after studying the exact Sangaku sources (985 files, incl. the full
docs/CAS.md). I built Lizard from source, loaded Sangaku, and ran the
relevant capabilities LIVE in the sandbox. Findings below are tested, with
one correction to the earlier plan and one capability gap to design around.

## Toolchain (verified working in-sandbox)
- Lizard builds from C (deps: `ds`, `gmp`); runs Sangaku via
  `cat src/prelude.lisp FILE.lisp | lizard`.
- Sangaku: proof-carrying CAS, ~21k lines Lisp / 239 modules. Ran examples
  139 (certificates) and 161 (linear algebra) successfully.

## THE keystone for dynsys: exact symbolic linear algebra (cas/linalg.lisp)
This is the highest-value module and I verified it live. Matrices are lists
of rational rows. Public API:
- `(mat-charpoly A)`            -> characteristic polynomial (Faddeev-LeVerrier)
- `(mat-eigenvalues A)`         -> EXACT eigenvalues (roots of charpoly): rationals, surds, or i
- `(mat-eigenvalues->string A)` -> human-readable eigenvalues
- `(mat-det A)`, `(mat-inverse A)`, `(mat-solve A b)`, `(trace M)`
- `(mat-charpoly->string A)`
- certificates: Cayley-Hamilton p(A)=0, A*A^-1=I, det cross-check,
  eigenvalues back-substituted, rational eigenvalues via det(A-lam I)=0.

LIVE results I obtained (exactly, not numerically):
- saddle [[1,0],[0,-1]]      -> eigenvalues -1, 1          charpoly x^2-1
- center [[0,-1],[1,0]]      -> eigenvalues +i, -i         charpoly x^2+1
- stable spiral [[-1,-1],[1,-1]] -> -1 +/- i               charpoly x^2+2x+2
- 3x3 diag(2,-3,5)           -> eigenvalues 2,5,-3          charpoly x^3-4x^2-11x+30

Why this matters: a center's eigenvalues are EXACTLY +/- i. Numerically they
come out as +/-1e-8 i and you cannot tell a true center / Hopf point from a
weak spiral. Sangaku certifies Re=0 exactly -- precisely what's needed to
nail a Hopf bifurcation and to classify equilibria without floating-point
ambiguity.

## Exact equilibria (verified)
- `(solve-poly p)` -> exact roots with multiplicity. x^2-2 -> +/-sqrt2 as
  algebraic numbers; x^2-3x+2 -> exact 1 and 2.
- `(count-real-roots-in p a b)` -> Sturm-exact count on an interval.
- Polynomial SYSTEMS: `(psys-solve F)` (Groebner) + `(psys-zero-dim? F nv)`
  + `(psys-eliminant F nv i)` -> exact multi-variable equilibria via
  elimination.
- bridge primitives: `(expr->poly e var)`, `(poly-deriv p)`, `(poly-eval p x)`.

## Certified differentiation (verified)
cas/diff-cert.lisp: `(diff e)` -> (value . certificate); `(derivative e)`;
`(certify e)` -> kernel-checked. Rules for +,*, sin, cos, exp, ln, neg,
recip, and the chain rule der_comp. A wrong derivative's certificate is
REJECTED by the dependent-type kernel (saw this live in example 139).

## Correction to the earlier plan + the one gap
EARLIER I implied Sangaku would hand us multivariate Jacobians directly.
The sources show single-variable `poly-deriv` and single-variable certified
`derivative`, the linear-algebra/eigenvalue pillar, and Groebner systems --
but NO ready-made multivariate partial-derivative or `jacobian` function.

Design consequence (clean, not a blocker): dynsys computes the JACOBIAN
ENTRIES itself -- it already has forward-mode AD (numerically exact at a
point) and its own expr IR -- then hands the resulting RATIONAL MATRIX to
Sangaku's mat-eigenvalues / mat-charpoly for EXACT eigenvalues + certified
classification. So:
  dynsys: vector field -> (its own) Jacobian at equilibrium, as rationals
  Sangaku: mat-eigenvalues(J) -> exact spectrum + certificate
This still delivers certified stability classification; only the matrix
assembly stays on the dynsys side (where it already is).

## Phase E features, re-prioritized by what the sources actually support
1. EXACT EIGENVALUES / certified stability at equilibria  [linalg.lisp]
   Highest value, lowest risk. Feed J (rationals) -> exact eigenvalues,
   classify (saddle/node/focus/center) with Re/Im known exactly, including
   the Re=0 Hopf-critical case. Show charpoly + eigenvalues + "kernel
   certified" in the UI.
2. EXACT EQUILIBRIA for polynomial fields  [solve-poly / psys-solve]
   1-D via solve-poly; 2-D+ via Groebner elimination. All roots, exact,
   verified -- strictly better than the numeric scan.
3. ANALYTIC NULLCLINES  [expr->poly + solve]
   Exact f_i=0 curves instead of marching-squares.
4. NORMAL-FORM / Hopf coefficients  [diff-cert + linalg]
   Build from exact derivatives + exact spectrum; turns numeric detection
   into certified classification.

## Bridge architecture (cas_bridge)
1. Serializer: tpcas AST/expr_ir -> Sangaku forms. Two targets: dense
   coefficient lists for univariate polynomials, and rational matrices
   (lists of rational rows) for Jacobians dynsys assembles.
2. Runner: write a generated .lisp, `cat prelude + script | lizard`, capture
   stdout, timeout + graceful fallback to the existing numeric path if the
   CAS/Lizard binary isn't present.
3. Parser: read back `(rat n)`, `(alg (min ...) (...))`, `+i`, coefficient
   lists, and proof-record s-exprs into dynsys numeric + display forms.

## Build/runtime notes (from getting it working here)
- Lizard needs ds (github:hydrastro/ds) and gmp; build with gmp as a SYSTEM
  include (-isystem) so its header doesn't trip Lizard's -Werror.
- The bridge shells out to a `lizard` binary; dynsys should locate it via
  $LIZARD or a configured path, and degrade gracefully when absent (the CAS
  features become "unavailable", the numeric tool keeps working).

## Licensing (unchanged, still the gating item for SHIPPING)
Sangaku's LICENSE currently grants no usage rights beyond reading until a
license is chosen. We can develop+test the bridge locally now; bundling or
shipping waits on the license you pick.

## Suggested first slice (one iteration, fully testable here)
Certified eigenvalues at equilibria:
- dynsys assembles J at a (numeric or exact) equilibrium as a rational matrix;
- serialize -> call mat-eigenvalues + mat-charpoly via Lizard;
- parse exact eigenvalues; classify with Re/Im known exactly;
- show in the UI next to the equilibrium, with the certificate flag.
Headless-testable: compare exact eigenvalues vs the current numeric eigen
path on known systems (saddle, center, stable focus, 3-D node).
