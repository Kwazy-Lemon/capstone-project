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

## Architecture

```text
Circuit Elements
       ↓
    NodeMap
       ↓
Element Stamping
       ↓
 MNA Assembler
       ↓
 Sparse MNA Matrix
       ↓
  Linear Solver
       ↓
 Circuit Results
```

## DC Benchmark

The current benchmark uses:

- Voltage source: 10 V
- R1: 1 kΩ
- R2: 2 kΩ

The unknown vector is:

```text
x = [Vin, Vout, I_Vs]^T
```

The assembled MNA system is:

```text
[ 0.001  -0.001   1    ] [Vin ]   [ 0 ]
[-0.001   0.0015  0    ] [Vout] = [ 0 ]
[ 1        0       0   ] [I_Vs]   [10 ]
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

## Testing

The simulator was tested at both the individual module level and the complete circuit level.

### 1. Module Import Test

All simulator modules can be imported successfully:

```text
All modules imported successfully.
```


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


### 3. Resistor Stamping Test

For a 1 kΩ resistor between `Vin` and `out`, the generated MNA contributions are:

```text
[(0, 0, 0.001),
 (1, 1, 0.001),
 (0, 1, -0.001),
 (1, 0, -0.001)]
```



### 4. Voltage Source Stamping Test

For a 10 V voltage source between `Vin` and `GND`, with branch-current index 2:

```text
([(0, 2, 1.0),
  (2, 0, 1.0)],
 (2, 10.0))
```


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

This confirms that the MNA system is assembled from the circuit parameters rather than being hard-coded for a single benchmark.

## Running

From the `circuit_simulator` directory:

```bash
python main.py
```

## Next Steps

- AC analysis
- Transient analysis using Backward Euler
- Nonlinear diode analysis using Newton iteration
- Additional circuit elements
