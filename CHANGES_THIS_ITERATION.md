# dynsys — IFS coefficient sliders are now editable (per-coefficient)

Unzip at repo root, then `make clean && make && make run`. Green
"dynsys NEW-UI <date>" label = you're on this build.

## The fix

I had disabled ALL the coefficient sliders for any IFS that had a parameter
or a non-constant coefficient — too blunt. Now editability is decided
**per coefficient**:

- A coefficient that's a plain constant (e.g. the fern's 0.85, -0.04, 1.6)
  gets a **live, editable slider** — drag it and the attractor updates.
- A coefficient that's a parameter expression (e.g. the spiral's
  `s*cos(theta)`) is shown as a **disabled slider in orange** that tracks
  its live value — you change those via the parameter sliders.

So a constant coefficient is editable **even when the system also has
parameters elsewhere**. A fully-constant IFS (fern, Sierpinski, dragon) is
now entirely editable, including the negative coefficients.

## A bug I caught before shipping

Negative literals like `-0.04` parse as `0 - 0.04` (a subtraction node), not
a bare number — so my first "is it a literal?" check wrongly flagged them
non-editable, which would have made the fern read-only. Fixed by testing the
real criterion: a coefficient is editable iff it references **no parameter**
(checked by walking the expression for any variable). Locked with a test
(`make test-ifslit`): fern coefficients incl. negatives = editable; the
spiral's `s*cos(theta)` = parameter-driven, its `e,f,p` = editable.

Editing one constant coefficient of a parametrized IFS keeps that edit while
the parameter-driven coefficients keep tracking their parameters.

## Verification
- All 4 C++ TUs compile with ZERO warnings (-O2 -Wall -Wextra).
- make test: 25 checks/suites pass (added **test-ifslit**).
