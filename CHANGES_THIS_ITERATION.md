# dynsys — guaranteed zoom controls + food-chain homoclinic investigation

Unzip at repo root, `make clean && make && make run`. (git re-initialized after a
sandbox reset; recovered intact from the shipped zip.)

## Phase-plane zoom: explicit buttons (addresses the recurring zoom-out report)
After auditing both the 2D and 3D zoom paths (the math is correct: wheel-down
gives an enlarging factor, manual bounds clamp only at +/-1e12, equal-aspect only
adjusts aspect, a left-click seed does NOT reset bounds), the most robust fix for
a mouse/trackpad-dependent issue is a guaranteed path that doesn't rely on the
wheel at all. The phase-plane control strip now has explicit
**[zoom out] [zoom in] [reset view]** buttons:
- zoom out enlarges the span 1.4x about the current view centre,
- zoom in shrinks it,
- reset view re-fits to the data (auto-bounds).
They set MANUAL bounds, so the auto-fit (which only ever grows the view) can no
longer fight a deliberate zoom-out. Hint text updated.

(Mouse-wheel zoom itself is unchanged and was verified correct by inspection;
the buttons are an additional, deterministic path. The interactive feel still
needs hands-on confirmation since the headless sandbox has no pointer.)

## Roadmap: find_homoclinic applied to the food chain
The food chain's cycle period diverges to ~300 inside the period-doubling bubble
(K~0.7-0.9), which suggested an accumulating homoclinic. Investigated with the
new find_homoclinic: it correctly reports this is NOT a simple homoclinic to an
interior saddle -- the reduced (x3=0) predator-prey subsystem has no positive
coexistence fixed point at these parameters, and the finder finds no same-saddle
excursion off the prey-only boundary equilibrium. The diverging period is instead
consistent with a near-BOUNDARY heteroclinic cycle (the orbit lingers near the
axial equilibria where a species is rare) -- a well-known tritrophic-food-chain
mechanism. This is the honest outcome: find_homoclinic distinguishes "no
same-saddle homoclinic here" from a false positive. (See docs/EXAMPLE_FOOD_CHAIN.md.)

## Verification
- All 4 C++ TUs compile ZERO warnings; the REAL GLFW/OpenGL binary builds.
- CAS green; find_homoclinic_smoke and validation_smoke pass. The analysis
  library is unchanged this iteration (UI + docs only), so the full suite is as
  previously verified.
