# Python Circuit Simulator

A modular Python implementation of a SPICE-like circuit simulator based on Modified Nodal Analysis (MNA).

## Current Scope

The current implementation supports:

- Node mapping
- Resistor modeling
- Independent voltage source modeling
- MNA matrix assembly
- Sparse matrix construction using COO and CSC formats
- DC circuit solving
- Analytical benchmark validation
- Reusable circuit parameter testing

AC analysis, transient analysis, nonlinear device models, and netlist parsing are not implemented in this version.

## Project Structure

```text
circuit_simulator/
├── README.md
├── main.py
├── node_map.py
├── resistor.py
├── voltage_source.py
└── mna_assembler.py
```

| File | Responsibility |
|---|---|
| `main.py` | Builds the example circuit, runs DC analysis, prints the solution, and checks the output voltage |
| `node_map.py` | Maps non-ground node names to matrix indices; ground has no voltage unknown |
| `resistor.py` | Defines a resistor and returns its conductance-matrix contributions |
| `voltage_source.py` | Defines an ideal voltage source and returns matrix entries and its voltage constraint |
| `mna_assembler.py` | Registers elements, assigns voltage-source current indices, and assembles the sparse MNA system |

## Architecture and Program Flow

The program follows this sequence:

1. `main.py` creates the circuit elements and registers them with the assembler.
2. `NodeMap` assigns indices to non-ground nodes. Ground is the reference node at 0 V and is not included in the unknown vector.
3. Each ideal voltage source receives an additional unknown for its branch current.
4. Each element's `stamp()` method returns the matrix entries contributed by that element.
5. `MNAAssembler` combines the contributions into a sparse matrix and a right-hand-side vector.
6. The example converts the small matrix to a dense array, solves the linear system, and compares the result with an analytical calculation.

```text
Circuit elements
      |
      v
  NodeMap ----------------> Node indices
      |
      v
Element stamping
      |
      v
 MNAAssembler ------------> Sparse matrix A and RHS b
      |
      v
 Solve A x = b
      |
      v
Node voltages and source currents
      |
      v
Analytical validation
```

## MNA Formulation

The assembled system has the form:

```text
A x = b
```

For the benchmark, the unknown vector is:

```text
x = [Vin, Vout, I_Vs]^T
```

It contains non-ground node voltages first, followed by currents through ideal voltage sources.

### Resistor contribution

For a resistor with resistance `R`, its conductance is `g = 1/R`. The resistor contributes:

- `+g` to the diagonal entry for each non-ground terminal
- `-g` to the off-diagonal entries when both terminals are non-ground
- No matrix entry for ground, because ground has no unknown index

### Ideal voltage-source contribution

An ideal voltage source imposes the constraint:

```text
V(+) - V(-) = source_voltage
```

Because the source current cannot be represented by a conductance alone, MNA adds a separate current unknown. The current appears in the connected nodes' KCL equations, and the voltage constraint adds a row and column to the matrix.

## DC Benchmark

The benchmark uses a 10 V source and a two-resistor voltage divider:

- Voltage source: 10 V between `Vin` and ground
- `R1`: 1 kΩ between `Vin` and `out`
- `R2`: 2 kΩ between `out` and ground

The expected output voltage is calculated independently:

```text
Vout = Vs * R2 / (R1 + R2)
     = 10 * 2000 / (1000 + 2000)
     = 6.666667 V
```

The assembled MNA system for the node order used by the example is:

```text
[ 0.001  -0.001   1    ] [Vin ]   [ 0  ]
[-0.001   0.0015  0    ] [Vout] = [ 0  ]
[ 1        0       0    ] [I_Vs]   [10 ]
```

The Python implementation produces:

```text
Vin  = 10.000000 V
Vout = 6.666667 V
I_Vs = -0.003333 A
```

Analytical result:

```text
Expected Vout = 6.666667 V
Error          = 0.000000e+00 V
```

![DC benchmark result](https://github.com/user-attachments/assets/98f00f19-5777-4cea-87c5-13a726f99e4f)

The negative source-current value follows the sign convention used for the MNA branch-current unknown; the output voltage agrees with the analytical divider equation.

## Dependencies and Running

Requirements:

- Python 3
- NumPy
- SciPy

Install the dependencies if necessary:

```bash
python -m pip install numpy scipy
```

From the `circuit_simulator` directory, run:

```bash
python main.py
```

The current example converts the sparse matrix to a dense array because the benchmark is very small. A larger circuit would generally benefit from an appropriate sparse linear solver.

## Testing

The following manual checks were recorded during development. The screenshots below are retained as evidence of the original checks.

### 1. Module Import Test

All simulator modules can be imported successfully:

```text
All modules imported successfully.
```

![Module import test](https://github.com/user-attachments/assets/01a407fc-c380-4c79-86e6-726579ebadd0)

### 2. NodeMap Test

The node mapping was verified using:

```text
Vin  → 0
out  → 1
GND  → None
```

The number of non-ground nodes is:

```text
2
```

![NodeMap test](https://github.com/user-attachments/assets/0925d74d-ab3e-4b83-ad6e-726579ebadd0)

### 3. Resistor Stamping Test

For a 1 kΩ resistor between `Vin` and `out`, the generated MNA contributions are:

```text
[(0, 0, 0.001),
 (1, 1, 0.001),
 (0, 1, -0.001),
 (1, 0, -0.001)]
```

![Resistor stamping test](https://github.com/user-attachments/assets/f35eb65f-352e-4771-9a30-df76b8ea091b)

### 4. Voltage Source Stamping Test

For a 10 V voltage source between `Vin` and `GND`, with branch-current index 2:

```text
([(0, 2, 1.0),
  (2, 0, 1.0)],
 (2, 10.0))
```

![Voltage source stamping test](https://github.com/user-attachments/assets/731fdecc-3613-45bc-96a3-40194f88780e)

### 5. DC Benchmark Test

For:

```text
Vs = 10 V
R1 = 1 kΩ
R2 = 2 kΩ
```

the simulator produces:

```text
Vin  = 10.000000 V
Vout = 6.666667 V
I_Vs = -0.003333 A
```

The analytical error is:

```text
0.000000e+00 V
```

![DC benchmark test](https://github.com/user-attachments/assets/c4c0701f-89a9-434b-a6e2-18c87e79d294)

### 6. Reusability Test

The circuit parameters were changed without modifying the MNA assembly implementation:

```text
Vs = 10 V
R1 = 2 kΩ
R2 = 2 kΩ
```

The simulator produces:

```text
Vin  = 10.0 V
Vout = 5.0 V
I_Vs = -0.0025 A
```

The expected result is:

```text
Vout = 5.0 V
```

![Reusability test](https://github.com/user-attachments/assets/bb6b8646-c9cf-400f-b3bd-01c94cbe2765)

This confirms that the MNA system is assembled from circuit parameters rather than hard-coded for a single benchmark.

## Limitations and Next Steps

Current limitations:

- The circuit is constructed directly in Python; SPICE netlist parsing is not implemented.
- Only resistors and independent ideal voltage sources are modeled.
- Input validation for invalid resistance values and malformed element definitions is not comprehensive.
- The manual checks are not yet organized into an automated test suite.

Possible next steps:

- AC analysis
- Transient analysis using Backward Euler
- Nonlinear diode analysis using Newton iteration
- Additional circuit elements
