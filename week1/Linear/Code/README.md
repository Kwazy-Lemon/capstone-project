# Linear Analysis Code

This folder contains the MATLAB code developed for the Week 1 linear circuit analysis prototype.

The code is organized around the three linear circuit analyses studied in this week:

- **DC Analysis**
- **AC / Frequency-Domain Analysis**
- **Transient / Time-Domain Analysis**

The purpose of these scripts is to manually construct and solve the corresponding MNA equations and to understand the numerical procedures required for the final Python-based circuit simulator.

## Code Structure

| File | Description |
|------|-------------|
| `dc_analysis.m` | Manual MNA formulation and solution for DC analysis |
| `ac_analysis.m` | MNA formulation and frequency-domain circuit analysis |
| `transient_analysis.m` | MNA formulation and time-domain circuit analysis |

## Implementation Approach

The MATLAB prototypes are intentionally kept simple and manually constructed.

The main focus is to understand:

- MNA matrix construction
- Circuit element contributions
- System equation formulation
- Numerical solution
- Analysis-specific procedures for DC, AC, and transient analysis

These scripts are prototypes for algorithm verification and are not intended to represent the final automated circuit simulator.

## Results

The corresponding simulation results and plots will be stored in:

```text
../figures/
