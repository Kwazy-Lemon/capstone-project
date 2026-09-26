# Linear Circuit Analysis

## 1. Objective

The objective of this section is to study the formulation and solution of linear circuit analysis using **Modified Nodal Analysis (MNA)**.

A simple linear circuit will be manually formulated and analyzed in MATLAB to understand the fundamental procedures required for the final Python-based circuit simulator.

The study will cover:

- DC Analysis
- AC / Frequency-Domain Analysis
- Transient / Time-Domain Analysis

The MATLAB prototype is intended to validate the underlying algorithms rather than serve as a complete or automated circuit simulator.

---

## 2. MNA Formulation

For a linear circuit, the circuit equations will be formulated using Modified Nodal Analysis.

The prototype will focus on understanding:

- Circuit node definition
- MNA unknown variables
- Matrix construction
- Right-hand-side vector construction
- Numerical solution of the resulting system

The contribution of individual circuit elements to the MNA system will be studied through a simple manually constructed circuit.

---

## 3. DC Analysis

The DC analysis will investigate the steady-state behavior of a linear circuit.

The main steps include:

1. Define the circuit and its nodes.
2. Formulate the MNA equations.
3. Construct the system matrix and right-hand-side vector.
4. Solve the resulting linear system.
5. Obtain the circuit voltages and currents.

The MATLAB implementation will be kept simple and manually constructed to focus on understanding the analysis procedure.

---

## 4. AC / Frequency-Domain Analysis

The AC analysis will extend the MNA formulation to the frequency domain.

The study will focus on:

- Frequency-dependent circuit behavior
- Representation of capacitors and inductors in the frequency domain
- Construction of the frequency-domain MNA system
- Solving the circuit response at different frequencies

The resulting frequency response will be used to examine the behavior of the circuit across the selected frequency range.

---

## 5. Transient / Time-Domain Analysis

The transient analysis will investigate the time-domain behavior of the circuit.

The main topics include:

- Time discretization
- Time-step progression
- Treatment of capacitors and inductors
- Use of previous time-step information
- Repeated solution of the circuit equations over time

The objective is to understand what additional steps are required for transient analysis when the relevant circuit matrices have already been established.

---

## 6. Expected Results

The linear prototype is expected to provide:

- A manually constructed MNA formulation
- A basic DC analysis
- A basic AC / frequency-domain analysis
- A basic transient / time-domain analysis
- A clearer understanding of the numerical procedures required for the final Python simulator

The results and observations from the MATLAB prototype will be used to guide the subsequent design of the Python implementation.
