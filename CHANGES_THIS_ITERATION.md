# dynsys — IFS maps are now inspectable & editable (the IFS "variables")

Unzip at repo root, then `make clean && make && make run`. Green
"dynsys NEW-UI <date>" label = you're on this build.

## What I'd missed, now fixed

You meant the IFS's own contents — its maps and their coefficients — should
be visible and editable, the way an ODE shows its equations and state
variables. An IFS was a black box: load it, see the attractor, but nothing
showed the maps. Now there's an **"IFS maps"** panel (Setup tab, shown when
an IFS is loaded) listing every affine map and its coefficients:

    x' = a x + b y + e      y' = c x + d y + f      (p = selection weight)

- For a plain (constant-coefficient) IFS — fern, Sierpinski, dragon — every
  coefficient is **editable**: drag a value (double-click to type) and the
  attractor updates live. You can also **+ add map / - remove last** to grow
  or prune the system.
- For a parametrized IFS (coefficients are expressions like `s*cos(theta)`),
  the panel shows the **evaluated values read-only** and tells you to adjust
  the parameters above — so a number-edit can't silently fight an
  expression. (Detected automatically at compile.)

Try it: load **"Barnsley fern (IFS)"**, open **IFS maps**, and drag map 2's
`d` from 0.85 downward — the fern collapses into a stunted plant in real
time. Load **"Spiral IFS (parametrized)"** and the panel is read-only with
the theta/s sliders driving it instead.

## Verified
- The editable path feeds edits straight to the chaos game: changing the
  fern's `d` (0.85 -> 0.5) shrinks its height from 10.0 to 3.2, and adding a
  map still yields a valid bounded attractor (checked headlessly + rendered).
- Editability is detected correctly: all-constant IFS = editable; any
  parameter expression = read-only.
- All 4 C++ TUs compile with ZERO warnings (-O2 -Wall -Wextra).
- make test: 24 checks/suites pass.

## Note on the two IFS knobs
- **IFS maps panel** = the maps themselves (this iteration).
- **Parameters panel** = parameters that *drive* expression coefficients
  (previous iteration). Both can apply to the same system; they're
  complementary.
