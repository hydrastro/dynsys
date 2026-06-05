# Should you build an FP / non-symbolic module? — recommendation

Short answer: **yes, but a small, sharply-scoped one, and most likely it
belongs in dynsys (or as a thin Lizard primitive), NOT inside Sangaku.**

## Why Sangaku should stay exact-only
Sangaku's whole identity is "every positive result carries a
machine-checkable certificate." Its roots come back as algebraic numbers
(even cplxroots/polyroots), there is deliberately no `->inexact`/decimal
layer, and that purity is the moat. Adding a floating-point evaluator inside
Sangaku would dilute the one thing that makes it special and create a
two-tier trust story inside a system whose pitch is uniform trust. Keep
Sangaku exact.

## Where an FP layer genuinely earns its place
There are three distinct jobs an FP/numeric layer could do. Only the first
is clearly worth it for dynsys today:

1. **Decimal REALIZATION of exact results (worth it, small).**
   Sangaku says an eigenvalue is `(1+sqrt 33)/2` or a root is a specific
   algebraic number. For a plot axis, a slider readout, or a "≈ 3.372"
   tooltip, you need a decimal. The exact value already pins the answer;
   you just want N correct digits. This is a *narrow* capability:
   "evaluate this exact algebraic number / rational / surd to k decimal
   digits, with a guaranteed error bound." Because the exact value is known,
   the FP result can even be *certified* (interval around the true value),
   so it doesn't break the Sangaku spirit — it's a presentation layer.

   Best home: a small helper, either
   - in **dynsys** (it already does all the floating-point plotting), taking
     the parsed exact form from the bridge and producing a double — simplest,
     and where the need actually lives; or
   - as a thin **Lizard** primitive `(alg->decimal a k)` if you want it
     reusable across Sangaku clients. Either is ~a day of work.

2. **Numerical root-finding / ODE integration (NOT needed — dynsys already
   has it).** dynsys already integrates ODEs (RK4/DOPRI/etc.), finds fixed
   points numerically, and computes numeric eigenvalues. Duplicating that in
   a new module is wasted effort. The division of labor is already right:
   dynsys = fast floating-point dynamics; Sangaku = exact analysis.

3. **Arbitrary-precision/interval arithmetic as a hedge (defer).** Useful
   only if you later want certified decimals to many digits or validated
   numerics for stiff cases. Lizard already links gmp, so a future
   `mpfr`-style layer is feasible — but it's speculative until a concrete
   need appears. Don't build it now.

## Concrete recommendation
- **Do** add a tiny "exact → decimal (k digits, with error bound)"
  realizer. Put it in **dynsys** first (closest to the need: axis labels,
  readouts, tooltips next to the exact eigenvalues/roots the bridge
  returns). Promote it to a Lizard primitive only if a second Sangaku client
  needs it.
- **Don't** put floating point inside Sangaku, and **don't** build a numeric
  ODE/rootfinder module — dynsys already owns that layer.
- **Defer** arbitrary-precision/interval arithmetic until a real use case
  (certified many-digit output, validated stiff numerics) shows up.

This keeps the architecture clean: Sangaku exact and certified, dynsys
floating-point and fast, and a thin realization seam between them that
turns exact answers into the decimals a GUI needs — without compromising
either side.
