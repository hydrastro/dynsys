# dynsys — IFS map coefficients are now real SLIDERS

Unzip at repo root, then `make clean && make && make run`. Green
"dynsys NEW-UI <date>" label = you're on this build.

## Fix

Last iteration I gave the IFS maps panel drag-to-scrub number fields, not
sliders — that's not what you wanted. Now every coefficient is a proper
**slider** (visible track + handle, click anywhere on the track), with a
sensible range per coefficient:
- a, b, c, d (linear / rotation part): slider over [-1, 1]
- e, f (translation): slider over [-5, 5]
- p (selection weight): slider over [0, 1]
(If a value sits outside its default window, that slider's range widens to
include it, so nothing is unreachable.)

Layout: two sliders per row at 150px each, grouped under "map 1", "map 2",
... so they're wide enough to actually use.

- **Constant-coefficient IFS** (fern, Sierpinski, dragon): the sliders edit
  the maps directly and the attractor updates live. + add map / - remove
  last still there.
- **Parametrized IFS** (coefficients are expressions): the coefficient
  sliders are shown disabled and track the live values as you drag the
  PARAMETER sliders above — so you can watch how, say, theta drives each
  a/b/c/d. (The editable knobs in that case are the parameters.)

## Try it
Load **"Barnsley fern (IFS)"**, open **IFS maps**, and slide map 2's `d`
down from 0.85 — the fern shrinks live. Load **"Spiral IFS
(parametrized)"** and slide `theta` in Parameters; the coefficient sliders
move to show the rotation feeding through.

## Verification
- All 4 C++ TUs compile with ZERO warnings (-O2 -Wall -Wextra).
- make test: 24 checks/suites pass (UI swap; the map/eval logic is the same
  validated code from last iteration).
