# dynsys — food-chain two-parameter study (codim-2 roadmap exercise)

Unzip at repo root, `make clean && make && make run`. (git re-initialized after a
sandbox reset; work recovered intact from the shipped zip.)

## Food chain now has a SECOND parameter (enables the two-parameter tools)
The "Food chain (3 species)" preset previously exposed only K. It now also
exposes **d** = predator death rate (default 0.4, range [0.1, 1.2]). With two
parameters the cycle-bifurcation CURVE tools become usable on this system:
seed a cycle and use "Trace period-doubling curve" / "Trace Neimark-Sacker
curve" to follow the locus in the (K, d) plane, with codim-2 detection
(fold-flip / PD-NS / cusp) along it. (No collision between the parameter `d` and
the dx1/dx2/dx3 derivative prefixes -- verified the preset still parses and
renders.)

## Two-parameter study results (verified headless)
- pd_curve traces the period-doubling locus in (K, d): e.g. from
  (K=0.70, d=0.42) to (K=0.64, d=0.44). The codim-2 detector runs cleanly with
  NO false positives on that segment.
- A 2-parameter sweep maps the period-doubling BUBBLE: it shifts to higher K as
  d increases (first-PD-K ~0.50 at d=0.30 -> ~0.92 at d=0.54) with width ~0.4;
  the two PD branches run roughly parallel over d in [0.3, 1.0], so the codim-2
  cusp where they would meet lies outside that window. Toward LOW d the cascade
  deepens (period 5 -> 6 -> 8 -> 10+) toward chaos.
- The cycle period is very long inside the bubble (~300 near K=0.7-0.9), dropping
  to ~55-80 for K>=1.1: the period-doublings accumulate near a HOMOCLINIC
  connection, which is why the period diverges there. This makes the food chain
  a genuine, non-textbook stress test of the long-period collocation/Floquet
  path (which traced it correctly).

See docs/EXAMPLE_FOOD_CHAIN.md for the full workflow and the bifurcation
inventory (Hopf + fold-of-cycles + a four-PD bubble).

## Verification
- All 4 C++ TUs compile ZERO warnings; the REAL GLFW/OpenGL binary builds and
  renders the 2-parameter preset.
- CAS green; codim2_cycle_smoke passes. The analysis library is unchanged this
  iteration (preset + docs only), so the full suite is as previously verified.
