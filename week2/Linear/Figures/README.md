# Week 2 — RC Transient Results

## Changes from Week 1

The initial sample is now stored correctly as Vout = 0 V at t = 0.
The simulation loop starts at the second sample, so the first calculated
voltage corresponds to t = 10 us.

The updated code also compares the numerical result with the exact RC
response and reports the maximum voltage error.

## Circuit Settings

- Input voltage: 10 V
- R1: 1000 ohms
- R2: 2000 ohms
- Capacitance: 100 nF
- Time step: 10 us
- Simulation duration: 5 ms
- Method: MNA with Backward Euler

## Expected Results

- Initial output voltage: 0 V
- Output voltage at 10 us: approximately 0.8696 V
- Steady-state output voltage: approximately 6.6667 V
- Time constant: approximately 66.67 us
- Maximum voltage error: approximately 0.1733 V

Backward Euler gives a slightly lower output than the exact solution
during the rise. A smaller time step reduces this error.

## Input and Output Voltage

The input stays at 10 V. The output rises from 0 V and approaches
6.6667 V. This plot shows the full 5 ms simulation.

<img width="966" height="807" alt="image" src="https://github.com/user-attachments/assets/eb50898b-b199-49a6-8fc1-b788cdbffb9a" />


## Comparison with the Exact Solution

This plot shows the first 300 us. It compares the Backward Euler
result with the exact RC response and makes the initial sample
and numerical error easier to see.

<img width="972" height="801" alt="image" src="https://github.com/user-attachments/assets/5a5a0e97-0b0e-4714-8b60-a70407c564fe" />

