# Roadmap: closing the remaining MatCont gaps

Empirical state (verified by grep + passing tests), Jun 2026.

## Already shipped & tested (NOT gaps — keep regression-guarded)
- Limit-cycle continuation: collocation BVP + pseudo-arclength + adaptive mesh +
  simulation self-seed.            [lpc_arclength_smoke, lc_colloc_smoke, lc_selfseed_smoke]
- Floquet multipliers via the variational monodromy; PD/NS *points* detected and
  bracketed on a branch.           [cycle_pdns_smoke]
- Codim-2 points on equilibrium curves: BT (a,b), GH (l2 sign), cusp (cubic c),
  zero-Hopf, Hopf-Hopf.            [bt_codim2_smoke, codim2_coeffs_smoke, zhhh_smoke]
- Homoclinic single-solve + continuation (tangent predictor). [homoclinic_smoke, homoclinic_cont_smoke]
- LPC (fold-of-cycles) curve in two parameters.  [lpc_arclength_smoke / lpccurve]
- Branch points + branch switching.  [branch_switch_smoke]

## GENUINE GAPS vs MatCont (the roadmap, in priority order)

### P1 — Two-parameter codim-1 CYCLE-bifurcation curves
MatCont follows these as curves in (p,q); we currently only detect the points on
a single branch. Each is: detect the point, then continue it by augmenting the
cycle BVP with the bifurcation condition and pseudo-arclength in (p,q).
  1. PD curve  (period-doubling locus): augment cycle BVP with the PD test
     (a real Floquet multiplier = -1) and continue in two params.
  2. NS curve  (Neimark-Sacker / torus locus): augment with the NS condition
     (a complex multiplier pair on the unit circle, |mu|=1).
  3. LPC curve already exists; verify/strengthen it sits in the same framework.

### P2 — Branch point of cycles (BPC) + cycle branch switching
Detect a singularity of the cycle bordered system (a genuine branch point of the
periodic-orbit family) and switch onto the crossing branch, mirroring the
equilibrium-branch switch_branch.

### P3 — Connecting-orbit breadth
  - Heteroclinic connections (saddle A -> saddle B) by a truncated BVP with
    projection BCs at both ends (generalizes solve_homoclinic).
  - Saddle-node homoclinic (HSN) / homoclinic to a non-hyperbolic equilibrium.

### P4 — Lin's method for the stiff Bogdanov-Takens homoclinic
Still experimental. Specific unsolved piece: a SIGNED Lin gap that brackets the
connection needs section geometry where both the forward (unstable) and backward
(stable) manifolds cross the SAME section cleanly; the returning loop crosses a
far-point section multiple times and the crossing selection is not yet robust.
solve_homoclinic / continue_homoclinic remain the validated path meanwhile.

## Method note
All P1 curves share one design: a bordered/augmented system
  [ cycle collocation BVP ;  bifurcation scalar condition ;  arclength ] = 0
solved by Newton with a pseudo-arclength predictor in the two active parameters.
We already have the cycle BVP (CycleCtx/cyc_newton), the Floquet machinery for
the conditions, and a working pseudo-arclength corrector for equilibria and
cycles — so P1 is mostly assembling existing pieces.

---

## Session refinement (mapping the 5 named areas to their OPEN sub-items)

The five headings (limit-cycle continuation, codim-1 curves, codim-2 points,
homoclinic continuation, Floquet/monodromy) are each largely DONE. The genuine
open work INSIDE them, in priority order:

A. [homoclinic] HETEROCLINIC connections — saddle A -> saddle B by a truncated
   BVP with projection BCs at BOTH ends. Direct generalization of the working
   solve_homoclinic. HIGH value, tractable now. ** START HERE. **
B. [limit cycle] BRANCH POINT OF CYCLES (BPC) — detect a singularity of the
   cycle bordered system along a branch and report it (switch later).
C. [codim-2] codim-2 points ON cycle-bifurcation curves — intersections like
   LPC-PD, PD-PD (period-doubling of the doubled cycle), and the strong
   resonances; detect where two cycle test functions vanish together.
D. [homoclinic] Lin's method for the stiff Bogdanov-Takens homoclinic — still
   experimental; the signed two-sided Lin gap on a shared section is the
   specific unsolved piece.

This session: implement A (heteroclinic connections) end to end with a
regression test, then proceed to B if budget allows.

---

## D (Lin's method) — diagnostic outcome

Investigated thoroughly. KEY FINDING: for the Bogdanov-Takens normal form
x'=y, y'=b1 + b2 y + x^2 + x y with b2 < 0 (the regime usually quoted for "the
BT homoclinic"), forward integration of the saddle's unstable manifold does NOT
return to that saddle — it connects to the OTHER equilibrium (a HETEROCLINIC,
which solve_heteroclinic now handles correctly) or settles on a nearby limit
cycle. This is exactly why a one-sided return test and a two-sided section gap
both fail here: there is no homoclinic-to-the-same-saddle to bracket in that
regime, and the backward stable manifold escapes to infinity rather than
crossing any section built from the unstable excursion (verified numerically).

OUTCOME (shipped): rather than a fully general Lin solver for the stiff case,
lin_homoclinic now DIAGNOSES what the unstable manifold actually does and
reports it honestly:
  - connects to a different equilibrium  -> "heteroclinic, use solve_heteroclinic"
  - settles on a cycle / escapes         -> "did not return; no homoclinic"
  - genuinely returns near the saddle with small gap -> located approximately.
It never emits a false positive. The validated connecting-orbit paths remain
solve_homoclinic / continue_homoclinic (homoclinic) and solve_heteroclinic
(saddle-to-saddle). A fully general Lin's method for the stiff same-saddle BT
loop (where a homoclinic and a limit cycle nearly coincide) remains the open
research item; it requires a genuine two-point BVP with the gap as an unknown,
not a shooting test function.
