# Week 2 — Resistor–Diode DC Analysis

## Changes from Week 1

The Newton solver now checks both the solution update and the
nonlinear circuit residual before reporting convergence.

Voltage and current use separate tolerances. Each diode-voltage
change is limited to 0.1 V, and unsafe numerical values stop the solver.

The code also saves convergence diagnostics and only plots a
validated operating point when all convergence checks pass.

## Circuit Settings

- Supply voltage: 5 V
- Resistance: 1000 ohms
- Diode saturation current: 1e-9 A
- Thermal voltage: 0.02585 V
- Ideality factor: 1
- Initial diode voltage: 0.4 V
- Maximum iterations: 100

## Convergence Checks

- Voltage tolerance: 1e-9 V
- Current tolerance: 1e-12 A

The scaled update measures the change between successive solutions.
The scaled residual measures how well the new solution satisfies
the original nonlinear circuit equations.

Both values must be at most 1 for convergence.

If the solver fails, it reports a warning and the last iterate
instead of claiming a validated operating point.

## Expected Results

For the default settings:

- Status: Converged
- Iterations: 4
- Diode voltage: approximately 0.3966 V
- Diode current: approximately 4.6034 mA
- Source current: approximately -4.6034 mA

The negative source current follows the MNA sign convention:
the voltage source supplies current to the circuit.

## Newton Convergence Diagnostics

This plot shows the scaled update and scaled nonlinear residual
at each iteration. Both must reach or fall below the threshold of 1.

<img width="1008" height="804" alt="image" src="https://github.com/user-attachments/assets/74281c59-a198-46cf-aed0-455ca7d3f7f8" />


## Diode DC Operating Point

This plot compares the diode I-V curve with the resistor load line.
Their intersection gives the DC operating point.

The marked point is shown only after convergence has been confirmed.

<img width="956" height="810" alt="image" src="https://github.com/user-attachments/assets/50d6b23b-66de-44ec-bedd-0b44628c9846" />

