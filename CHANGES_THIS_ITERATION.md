# dynsys — root-cause fixes + 19 new presets

Unzip at repo root, then `make clean && make && make run`. Green
"dynsys NEW-UI" label = you're on this build.

I found the real cause behind several recurring complaints this time.

## 1. Logistic bifurcation "horizontal lines" — ROOT CAUSE FOUND & FIXED

The Logistic and Tent presets were secretly "embedded in 3D"
(x_next, y_next=x, z_next=z, i.e. 3 state variables). That single bad
choice caused THREE of your bugs:
- the bifurcation observed a muddled variable -> not the clean tree;
- the 3D bridge (which needs a 1D map) was permanently disabled;
- the 1D cobweb view wasn't offered.

Both are now GENUINE 1D maps (state x only). Result:
- Logistic bifurcation now shows the real period-doubling tree (verified:
  the map gives 1,2,4,8,... attractor values as r increases).
- The 3D bridge works again for the Logistic/Tent (it needs a 1D map).
- The 1D cobweb view is offered for them.

## 2. The 3D bridge being "always deactivated" — fixed by the above

It is valid for 1D maps; with the Logistic/Tent now genuinely 1D, select
one of them and the "3D bridge" toolbar button is enabled.

## 3. Bifurcation "Run" button overflowing — replaced

The floating in-view panel that overflowed is gone. The bifurcation
controls (parameter picker, from/to range, Run) now live in the TOP
TOOLBAR when the bifurcation view is active — compact and clip-free. The
empty-view text now points you to the toolbar button and to loading the
Logistic map.

## 4. Basin / fractal display bugs — hardening

- Fixed an asymmetric bounds guard in the fractal parameter-space compute
  (one of the two parameter writes wasn't size-checked).
- (From last round, still in place) grid computes force fixed-step RK4 and
  are debounced so panning/zooming doesn't freeze.
Note: the fractal view requires a 2D MAP — use "Complex quadratic"
(Mandelbrot/Julia) or "Gingerbreadman". The now-1D Logistic correctly does
not offer the fractal view.

## 5. NEW PRESETS (19 added; 39 total)

Phase-plane: Duffing, FitzHugh-Nagumo, Brusselator, Sel'kov glycolysis,
SIR epidemic, Rosenzweig-MacArthur predator-prey (limit cycle), undamped
pendulum.
1D maps (great for bifurcation/cobweb): Sine, Gauss/mouse, Cubic.
2D maps: Gingerbreadman, Chirikov standard map.
3D: Chua's circuit (double-scroll), Sprott B, Nose-Hoover.
All 19 were checked through the real expression engine (parse + lower +
iterate); every one parses and runs.

## Verification
- All 4 C++ TUs compile with ZERO warnings (-O2) against the real
  ImGui/GLFW/cglm/tpcas headers (integer ImTextureID, matching nix).
- make test: 12 suites all pass.
- New presets validated by parsing+running each through the IR.

## What I still can't verify (need your eyes)
I can't see the GUI, so for any remaining "bug out" in basins/fractals: a
screenshot (with the green NEW-UI label) plus which preset + view triggers
it would let me fix the exact pixels. If it hard-crashes, running
`gdb ./build/dynsys` and pasting the backtrace pinpoints it instantly.
