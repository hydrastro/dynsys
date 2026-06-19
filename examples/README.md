# dynsys examples

Each `.dyn` file describes one dynamical system in the small text format
documented in the top-level `README.md`. Load one in the GUI from the preset
list, or run it headless:

```
./build/dynsys --headless examples/lorenz.dyn --steps 10000
```

| File | System | Kind |
|---|---|---|
| `lorenz.dyn` | Lorenz attractor | ODE, 3-D chaotic |
| `rossler.dyn` | Rössler attractor | ODE, 3-D chaotic |
| `thomas.dyn` | Thomas' cyclically symmetric attractor | ODE, 3-D |
| `four_dimensional_demo.dyn` | A 4-D demo system | ODE, 4-D |
| `van_der_pol.dyn` | Van der Pol oscillator | ODE, limit cycle |
| `damped_pendulum.dyn` | Damped pendulum | ODE, 2-D |
| `lotka_volterra.dyn` | Lotka-Volterra predator-prey | ODE, 2-D |
| `saddle_separatrix.dyn` | A saddle with separatrices | ODE, 2-D |
| `henon.dyn` | Hénon map | discrete map |
