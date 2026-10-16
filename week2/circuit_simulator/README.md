# Python Circuit Simulator

A modular Python implementation of a small SPICE-like DC circuit simulator. It uses **Modified Nodal Analysis (MNA)** to convert a circuit into a linear system and uses NumPy/SciPy to represent and solve that system.

## Current Scope

This version implements:

- Mapping non-ground circuit nodes to matrix indices
- Ideal resistor modeling
- Independent ideal voltage source modeling
- MNA matrix and right-hand-side assembly
- Sparse matrix construction with SciPy
- A small DC benchmark and an analytical voltage-divider check

AC analysis, transient analysis, nonlinear device models, and netlist parsing are not implemented in this version.

## Project Structure

| File | Responsibility |
|---|---|
| `main.py` | Builds and solves the example circuit, prints results, and checks the output voltage against an analytical calculation |
| `node_map.py` | Maps node labels to matrix indices; ground is the reference node and has no unknown |
| `resistor.py` | Defines a resistor and returns its conductance-matrix contributions |
| `voltage_source.py` | Defines an ideal voltage source and returns its matrix and voltage-constraint contributions |
| `mna_assembler.py` | Collects elements, assigns voltage-source current indices, and assembles the sparse MNA system |

## How the Code Works

The program follows this sequence:

1. `main.py` creates the circuit elements.
2. `MNAAssembler.add_element()` registers each element and its terminal nodes.
3. `NodeMap` assigns indices to non-ground nodes. Ground is the zero-volt reference and is not included in the unknown vector.
4. `MNAAssembler.assign_branch_indices()` assigns an additional unknown to the current through each ideal voltage source.
5. Each element's `stamp()` method returns the matrix entries contributed by that element.
6. `MNAAssembler.assemble()` combines the entries into a SciPy sparse matrix and creates the right-hand-side vector.
7. `main.py` solves the small benchmark system and compares the output voltage with the voltage-divider formula.

## MNA Formulation

The assembled system has the form:

```text
A x = b
```

The unknown vector contains non-ground node voltages first, followed by ideal voltage-source currents. For the example circuit:

```text
x = [Vin, Vout, I_Vs]^T
```

A resistor with conductance `g = 1/R` contributes:

- `+g` to the diagonal entry for each non-ground terminal
- `-g` to the off-diagonal entries when both terminals are non-ground
- No entry for a ground terminal, because ground has no matrix index

An ideal voltage source adds a current unknown and the constraint:

```text
V(+) - V(-) = source_voltage
```

The source current appears in the connected nodes' KCL equations. The voltage constraint adds a row and column to the MNA system.

## DC Benchmark

The example is a resistor divider:

- Voltage source: 10 V between `Vin` and ground
- `R1`: 1 kΩ between `Vin` and `out`
- `R2`: 2 kΩ between `out` and ground

The expected output voltage is calculated independently:

```text
Vout = Vs * R2 / (R1 + R2)
     = 10 * 2000 / (1000 + 2000)
     = 6.666667 V
```

The assembled system for this node order is:

```text
[ 0.001  -0.001   1    ] [Vin ]   [ 0  ]
[-0.001   0.0015  0    ] [Vout] = [ 0  ]
[ 1        0       0    ] [I_Vs]   [10 ]
```

Expected output from `main.py`:

```text
Vin  = 10.000000 V
Vout = 6.666667 V
I_Vs = -0.003333 A
Expected Vout = 6.666667 V
Error         = approximately 0 V
```

The negative source-current result follows the branch-current sign convention used by the MNA equations; it does not indicate an incorrect output voltage.

## Dependencies

- Python 3
- NumPy
- SciPy

Install the required libraries if needed:

```bash
python -m pip install numpy scipy
```

## Running the Example

Run from the `circuit_simulator` directory:

```bash
python main.py
```

The example converts the sparse matrix to a dense array for the NumPy solve because the benchmark is very small. For larger circuits, a suitable sparse linear solver should be considered.

## Testing and Limitations

The existing project notes record module-import checks, node-map checks, resistor and voltage-source stamping checks, the DC benchmark, and a second voltage-divider parameter set. These are useful manual checks; they are not yet organized as an automated test suite.

Current limitations:

- The program constructs the circuit directly in Python; it does not parse SPICE netlists.
- Only resistors and independent ideal voltage sources are modeled.
- Input validation for invalid resistance values and malformed element definitions is not yet comprehensive.
- AC, transient, and nonlinear analyses remain future work.
