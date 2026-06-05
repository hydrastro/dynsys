# What's missing in Sangaku (from a dynamical-systems standpoint)

Audited against the exact sources (985 files incl. docs/CAS.md). Sangaku's
CAS is very complete on algebra, integration, summation, and exact linear
algebra. The gaps below are the ones that matter for a dynamical-systems
tool; I implemented the first group (see cas/dynsys.lisp).

## Implemented this round (were missing; now added in cas/dynsys.lisp)
These had NO equivalent anywhere in src/ before:
- `mpoly-deriv p i` — partial derivative of a multivariate polynomial w.r.t.
  variable i. (Sangaku had single-variable poly-deriv and certified
  expression `derivative`, but no multivariate partial.)
- `vf-jacobian F nv` / `jacobian-at F nv pt` — the Jacobian of a polynomial
  vector field, symbolic and evaluated-at-a-point (rational matrix).
- `equilibrium?`, `equilibrium-eigenvalues/-charpoly/->string`,
  `equilibrium-trace/-det` — exact linearization spectrum at an equilibrium,
  composing the new Jacobian with linalg.lisp's mat-eigenvalues.
- `vf-divergence` / `divergence-at` — phase-volume contraction (trace field).
- `mpoly-gradient`, `mpoly-deriv2`, `mpoly-hessian`, `hessian-at` — gradient
  and second-derivative/Hessian for normal-form and extremum work.
Verified on Lorenz: Jacobian at origin exactly [[-10,10,0],[28,-1,0],
[0,0,-8/3]]; eigenvalues -8/3 and (-11 +/- sqrt 1201)/2; divergence the
constant -41/3. Golden test cas_dynsys + example 392 added.

## Still missing / candidate future work (NOT done)
Ordered by value to a dynamical-systems CAS:

1. **Symbolic partial derivatives of NON-polynomial fields.** diff-cert
   differentiates elementary expressions in ONE variable with a kernel
   certificate, but there is no multivariate certified partial for fields
   with sin/exp/etc. (e.g. a pendulum's sin(theta)). The polynomial Jacobian
   I added covers polynomial/rational fields; transcendental fields would
   need a multivariate extension of diff-cert. Medium effort, high value.

2. **EigenVECTORS (not just eigenvalues).** linalg.lisp gives exact
   eigenvalues and the kernel (nullspace) via normalform.lisp, so exact
   eigenvectors are assemblable as ker(A - lambda I), but there is no
   packaged `mat-eigenvectors`. Needed for invariant manifolds / mode
   shapes. Low-medium effort (compose existing pieces).

3. **A real-part / half-plane decision for algebraic eigenvalues.** Sangaku
   returns complex eigenvalues as algebraic numbers but exposes no "sign of
   the real part" or Hurwitz/Routh count. (dynsys does exact Routh-Hurwitz
   on the charpoly on its own side; folding a certified version into Sangaku
   would let the CAS itself emit a stability verdict.) Medium effort.

4. **Resultant-based elimination helper for equilibria.** psys-solve
   (Groebner) handles polynomial systems, and psys-eliminant extracts a
   univariate eliminant, but a turnkey "all real equilibria of this 2-D/3-D
   polynomial field, exactly" wrapper would make exact equilibria one call.
   Low effort (wrap psys + solve-poly).

5. **Interval / certified-decimal realization.** By design Sangaku is
   exact-only (no ->inexact layer). Turning an exact algebraic eigenvalue
   into "k correct decimal digits with an error bound" is the FP-realizer
   discussed separately; best placed in dynsys (or a thin Lizard primitive),
   NOT in Sangaku, to preserve its exact-and-certified purity.

## Non-gaps (deliberately out of scope; don't add)
- Floating-point linear algebra / numeric ODE integration / numeric
  rootfinding: dynsys already owns the fast numeric layer; duplicating it in
  Sangaku would only blur the clean exact-vs-numeric split.
- The HoTT/cubical limitations in docs/LIMITATIONS.md are Lizard-kernel
  concerns, not CAS gaps, and are irrelevant to dynsys.
