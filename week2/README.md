# Week 2 — Prototype Refinement and Python Simulator Development

## 1. Week 2 Objectives

Week 2 builds on the MATLAB MNA studies and preliminary software design from Week 1. The main objectives are to improve numerical validation, implement the first modular Python circuit framework, and develop reusable solvers for DC, AC, and transient analysis.

The work is organized into four areas:

- Refining the MATLAB linear RC transient prototype.
- Improving convergence checks for the MATLAB resistor–diode DC prototype.
- Implementing node mapping, element stamping, and MNA assembly in Python.
- Developing and validating Python numerical solvers using analytical benchmarks.

---

## 2. Work Completed

### 2.1 Linear Circuit Analysis in MATLAB

The RC transient prototype uses **MNA with Backward Euler**. The initial output voltage is now stored as **0 V at t = 0**, and the simulation starts calculating new states at the second sample. This corrects the alignment between the time vector and the computed response.

The updated script also compares the numerical waveform with the exact RC response, reports the maximum voltage error, and exports plots and a CSV result table.

- Circuit: 10 V source, R1 = 1 kΩ, R2 = 2 kΩ, and C = 100 nF.
- Simulation: 10 µs time step over 5 ms.
- Expected steady-state output: approximately 6.6667 V.
- Time constant: approximately 66.67 µs.
- Maximum voltage error: approximately 0.1733 V for the selected time step.

See [Linear transient results](Linear/Figures/README.md) and [the MATLAB script](Linear/transient_analysis.m).

### 2.2 Nonlinear Diode DC Analysis in MATLAB

The resistor–diode prototype continues to use **MNA and Newton iteration**, with stronger convergence and numerical checks:

- Both the solution update and the original nonlinear circuit residual must satisfy their tolerances.
- Voltage and current use separate tolerances: 1e-9 V and 1e-12 A.
- Each diode-voltage change is limited to 0.1 V per iteration.
- Non-finite values and unsafe exponential arguments stop the solver.
- Convergence history is exported, and the operating-point plot is generated only after convergence is confirmed.

For the default 5 V source and 1 kΩ resistor, the documented result converges in **4 iterations**, with a diode voltage of approximately **0.3966 V** and a diode current of approximately **4.6034 mA**. The source current is approximately **−4.6034 mA**, following the MNA reference direction.

See [Nonlinear analysis details](Non_linear/README.md) and [the MATLAB script](Non_linear/SimpleDiodeCircuit.m).

### 2.3 Modular Python Circuit Framework

The first Python framework separates circuit modeling from MNA matrix assembly:

| Module | Responsibility |
|---|---|
| `node_map.py` | Maps node names to matrix indices and handles ground nodes. |
| `resistor.py` | Generates resistor conductance stamps. |
| `voltage_source.py` | Generates voltage-source stamps and assigns branch-current unknowns. |
| `mna_assembler.py` | Collects element contributions and builds a COO matrix, then converts it to CSC format. |
| `main.py` | Constructs a resistor-divider benchmark, solves it, and compares the output with the analytical result. |

The benchmark uses a **10 V source, R1 = 1 kΩ, and R2 = 2 kΩ**. It produces:

```text
Vin  = 10.000000 V
Vout = 6.666667 V
I_Vs = -0.003333 A
```

The matrix assembler produces a sparse matrix. The current demonstration converts it to a dense array and uses `numpy.linalg.solve`; connecting it to the sparse numerical solver is a next step.

See [Python circuit framework](circuit_simulator/README.md).

### 2.4 Python Numerical Solvers and Benchmarks

The `numerical-solver` directory provides a separate reference implementation for linear circuit analysis:

- **DC:** solves the static MNA system using a sparse linear solver.
- **AC:** solves the complex frequency-domain system at each frequency.
- **Transient:** uses fixed-step Backward Euler and reuses the matrix factorization across time steps.

The benchmark script checks two DC circuits, an RC low-pass frequency response, RC charging with five time steps, and an RL circuit for inductor sign conventions. It exports voltages, currents, equation residuals, plots, error tables, and a validation summary.

The circuits in this package use explicitly defined matrices. These solvers have not yet been connected to the element-based assembly in `circuit_simulator`.

See [Numerical solver documentation](numerical-solver/ReadMe.md).

---

## 3. Numerical Validation Results

The committed [validation summary](numerical-solver/figure_results/validation_summary.json) reports:

| Check | Saved result |
|---|---|
| Two DC analytical benchmarks | Passed |
| RC AC maximum complex error | Approximately 2.22e-16 |
| RC cutoff frequency | Approximately 159.15 Hz |
| Gain and phase at cutoff | Approximately −3.0103 dB and −45° |
| RC transient checks | Passed |
| Maximum transient error at 12.5 µs | Approximately 0.01144 V |
| Observed transient convergence order | Approximately 0.993, consistent with first-order Backward Euler |
| RL DC, AC, and transient checks | Passed |

The RC numerical-solver benchmark uses **R = 1 kΩ and C = 1 µF**, with a **5 V step** for transient analysis. Its parameters differ from the MATLAB RC prototype above.

The [time-step error table](numerical-solver/figure_results/timestep_error.csv) shows the maximum transient error decreasing from approximately **0.16999 V at 200 µs** to **0.01144 V at 12.5 µs**.

A MATLAB reference implementation is included for independent comparison. The saved summary marks the MATLAB comparison as **pending**; the Python analytical checks and the Python–MATLAB comparison are separate validation steps.

---

## 4. Directory Structure

```text
week2/
├── README.md
├── Linear/
│   ├── transient_analysis.m
│   └── Figures/README.md
├── Non_linear/
│   ├── SimpleDiodeCircuit.m
│   └── README.md
├── circuit_simulator/
│   ├── README.md
│   ├── main.py
│   ├── node_map.py
│   ├── resistor.py
│   ├── voltage_source.py
│   └── mna_assembler.py
└── numerical-solver/
    ├── ReadMe.md
    ├── linear_solvers.py
    ├── run_benchmarks.py
    ├── matlab_reference.m
    └── figure_results/
```

---

## 5. How to Run

### Python

Install the required libraries, then run the two Python entry points from the **repository root**:

```bash
python -m pip install numpy scipy matplotlib
python week2/circuit_simulator/main.py
python week2/numerical-solver/run_benchmarks.py
```

The first script prints the assembled matrix and DC benchmark result. The second runs the numerical checks and saves CSV files, PNG plots, and `validation_summary.json` in `week2/numerical-solver/results_python/`.

The existing `figure_results/` folder contains committed benchmark evidence. New benchmark runs write to `results_python/`.

### MATLAB

Set MATLAB's Current Folder to the corresponding directory and run:

| Current Folder | Command |
|---|---|
| `week2/Linear` | `transient_analysis` |
| `week2/Non_linear` | `SimpleDiodeCircuit` |
| `week2/numerical-solver` | `matlab_reference` |

Run the Python numerical benchmarks before `matlab_reference` so MATLAB can read the newly generated files in `results_python/`. MATLAB reference results and comparison tables are written to `results_matlab/`.

---

## 6. Current Status

Week 2 establishes an element-based Python DC framework and independently validated linear solver routines, while improving the MATLAB transient and nonlinear prototypes.

