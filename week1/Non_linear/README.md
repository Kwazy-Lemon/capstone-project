# R + Diode DC Analysis (MNA in MATLAB)

Solves the DC operating point of a resistor in series with a diode
using Modified Nodal Analysis (MNA) and Newton–Raphson iteration.

## Circuit

Simple circuit with voltage source Vs with a resistor in series with a Diode.

## Files
- `SimpleDiodeCircuit.m` — main script

## How to run
1. Open MATLAB.
2. Run `SimpleDiodeCircuit.m`.
3. Read the converged node voltages and currents in the console.

## Parameters (edit at top of script)
- `Vs` — source voltage [V]
- `R`  — resistance [Ω]
- `Is` — diode saturation current [A]
- `Vt` — thermal voltage [V]

## Method
- Unknowns: `[v1; v2; i_vs]`
- Diode linearized each iteration with conductance `gd` and current source `Ieq`.
- NR stops when update < `1e-10` or 100 iterations.

## Notes
- Exponent is capped at 40 to avoid overflow.
- Initial guess `x = [0;0;0]` works for forward bias.
- Sign of `i_vs` depends on chosen convention.
## DC operating point
<img width="1757" height="1158" alt="DC_operating point" src="https://github.com/user-attachments/assets/ff0ff31d-84cd-4ff8-a543-171a65c94346" />

## Newton Convergence
<img width="1787" height="1154" alt="Newton_convergence" src="https://github.com/user-attachments/assets/208ee95c-87ca-418e-9d48-3144d609f886" />

