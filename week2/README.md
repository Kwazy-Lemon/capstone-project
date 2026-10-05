# Week 2 – Python-Based Circuit Simulator

## Overview

This week focuses on transitioning from the MATLAB-based circuit analysis prototypes developed in Week 1 to a Python-based circuit simulation framework.

The main goal is to establish the basic software architecture and numerical foundation required for a reusable SPICE-like circuit simulator.

---

## Objectives

The main objectives for Week 2 are:

1. Develop a basic Python framework for circuit simulation.
2. Implement node management and circuit element interfaces.
3. Formulate and assemble the Modified Nodal Analysis (MNA) system.
4. Support sparse matrix representation for the circuit equations.
5. Validate the Python implementation against the MATLAB results developed in Week 1.
6. Prepare the framework for future DC, AC, transient, and nonlinear analyses.

---

## Current Progress

### 1. MATLAB Linear Circuit Model

A MATLAB-based linear circuit model was developed during Week 1 to verify the basic circuit formulation and numerical approach.

The MATLAB implementation serves as a reference for validating the Python implementation.

### 2. MATLAB Nonlinear Circuit Model

A nonlinear circuit model was also developed to investigate the treatment of nonlinear circuit elements.

The results will be used as a reference when nonlinear device models are integrated into the Python simulator.

### 3. Python Circuit Simulator

The Python implementation is being developed as a modular and reusable framework.

The planned architecture includes:

```text
Circuit
│
├── NodeMap
│
├── Elements
│   ├── Resistor
│   ├── Capacitor
│   ├── Inductor
│   └── Nonlinear Elements
│
├── MNA Assembly
│
└── Solver
