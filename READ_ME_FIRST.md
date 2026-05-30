# READ THIS FIRST

## Build it so your changes actually take effect
From the repo root, AFTER unzipping this over it:
```
make clean      # deletes stale build/ — this step is essential
make
make run        # or ./build/dynsys
```
An unzip can leave old object files newer than the new source, so without
`make clean` the build may silently reuse the old binary.

## Confirm you're on the new build
At the TOP of the window you should see a green label:
    dynsys NEW-UI  (full-window plot)
If you don't see it, you're running an old binary — `make clean && make`.

## To see nullclines + their direction arrows
Nullclines and auto fixed points exist only for ODEs (a map has no
continuous vector field). Load **Van der Pol** or **Lotka-Volterra** from
the preset list — NOT the logistic/Hénon maps — and you'll see:
- red/green nullclines with little flow-direction arrows along them,
- auto-found equilibria, classified, with stable/unstable manifolds,
- the speed-colored vector field.
Toggle any of these in the "Phase plane" section of the left panel.

## Quick tour
- The plot fills the whole window; controls are the top toolbar + the
  left "controls" panel (toggle with "<< panel").
- 2D: left-click adds an orbit through that point (orbits accumulate;
  "Clear orbits" removes them), left-drag pans, wheel or +/- zooms,
  double-click auto-fits. "Add custom orbit" lets you type exact ICs.
- 3D (only offered for genuinely 3D systems): drag rotates, wheel/+/- zoom.
