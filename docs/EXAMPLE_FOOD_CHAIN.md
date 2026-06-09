# Example: tritrophic food chain (Rosenzweig-MacArthur)

Load the preset **"Food chain (3 species)"** (Continuous 3D category).

## The model
Prey x1, predator x2, top predator x3, with Holling type-II functional responses:

    dx1 = 0.11*x1*x3 - x1*x2/(x1+0.17) + (1-x1)*x1
    dx2 = K*x1*x2/(x1+0.17) - 0.4*(x2*x3/(x2+0.42) + x2)
    dx3 = 0.11*x2*x3/(x2+0.42) - 0.008*x3

The single parameter **K** (predator conversion efficiency) sweeps 0 -> 5.

## What to try (exercises the new tools)
- **Equilibria + stability**: at the default K=2.5 the coexistence state is an
  unstable focus (the phase-plane view labels it and draws the x1'/x2' nullclines
  and the two axial saddles).
- **Limit-cycle continuation in K**: from a simulated cycle, continue in K. The
  branch is non-trivial -- it turns around a FOLD OF CYCLES (LPC), which the
  pseudo-arclength continuation handles.
- **Floquet / period-doubling**: with Floquet enabled, a PERIOD-DOUBLING is
  detected near K ~ 1.44 (a real Floquet multiplier crosses -1); beyond it the
  cycle period-doubles toward chaos. Non-trivial multipliers grow large
  (tens to hundreds), reflecting strongly unstable cycles on the route to chaos.
- **Two-parameter curves**: the preset exposes a SECOND parameter **d**
  (predator death rate, default 0.4, range [0.1, 1.2]). Seed a cycle (e.g. at
  K=0.7, d=0.4) and use "Trace period-doubling curve" to follow the PD locus in
  the (K, d) plane; the codim-2 detector flags any fold-flip / PD-NS / cusp
  points along it.

## Two-parameter findings (verified)
- pd_curve traces the PD locus in (K, d) -- e.g. from (K=0.70, d=0.42) to
  (K=0.64, d=0.44) -- and the codim-2 detector runs cleanly (no false positives)
  on that segment.
- A 2-parameter sweep shows the PD BUBBLE shifts to higher K as d rises
  (first-PD-K ~0.50 at d=0.30 -> ~0.92 at d=0.54) with width staying ~0.4; the
  two PD branches run roughly parallel over [0.3, 1.0], so the codim-2 cusp is
  outside that window. Toward LOW d the cascade deepens (period 5->6->8->10+)
  into chaos.
- The cycle period is very long inside the bubble (~300 near K=0.7-0.9) and
  drops to ~55-80 for K>=1.1: the period-doublings accumulate near a HOMOCLINIC
  connection (hence the diverging period) -- a genuine, non-textbook stress test
  of the long-period collocation/Floquet path.

## Notes
Verified headless: the preset loads, renders, and the limit-cycle continuation
reports `traced 64 cycles (turned around a fold of cycles)` with a period-
doubling located at K=1.4392 and the seed cycle period ~54.

## Bifurcation inventory (verified, K in [0,5])
Mapped by limit-cycle continuation + Floquet and by an orbit-diagram scan:
- **Hopf**: the coexistence equilibrium is an unstable focus across this range
  (a Hopf has already occurred as K rose into [0,5]); a stable cycle exists.
- **Fold of cycles (LPC)**: the limit-cycle branch turns around a fold during
  pseudo-arclength continuation in K.
- **A period-doubling BUBBLE** in K ~ [0.67, 1.13] -- the cycle period-doubles
  and then un-doubles as K sweeps through:
    * PD1  (period-1 -> period-2)  near K ~ 0.68
    * PD2  (period-2 -> period-4)  near K ~ 0.76
    * a period-4 / weakly-chaotic band around K ~ 0.85-1.0
    * reverse PD2 (period-4 -> period-2) near K ~ 1.04
    * reverse PD1 (period-2 -> period-1) near K ~ 1.12
  i.e. at least FOUR period-doubling bifurcations bound the bubble.
So the system shows (at least): one Hopf, one fold-of-cycles, and a four-PD
bubble -- a good multi-bifurcation showcase. Higher-K and two-parameter sweeps
likely reveal more (Neimark-Sacker, branch points, codim-2 fold-flip / PD-NS).
