"""Run a small DC benchmark for the MNA circuit simulator.

The example is a 10 V source driving a two-resistor voltage divider. It shows
how circuit elements are registered, how the MNA system is assembled, and how
the numerical result is checked against the voltage-divider equation.
"""

import numpy as np

from resistor import Resistor
from voltage_source import VoltageSource
from mna_assembler import MNAAssembler


def main():
    """Build, solve, and validate the example DC circuit."""
    # Create an assembler to collect circuit elements and map their nodes.
    assembler = MNAAssembler()

    # Define a 10 V source and a 1 kΩ / 2 kΩ resistor divider.
    R1 = Resistor("R1", "Vin", "out", 1000.0)
    R2 = Resistor("R2", "out", "GND", 2000.0)
    Vs = VoltageSource("Vs", "Vin", "GND", 10.0)

    # Register elements before assembly so all node indices can be assigned.
    assembler.add_element(R1)
    assembler.add_element(R2)
    assembler.add_element(Vs)

    # Assemble A and b for the linear system A @ x = b.
    # x contains node voltages followed by voltage-source branch currents.
    A, b = assembler.assemble()

    # This benchmark is small, so convert the sparse matrix to dense form for
    # NumPy's solver and readable output. Larger systems should use a sparse
    # linear solver to avoid the cost of dense storage and factorization.
    A_dense = A.toarray()
    x = np.linalg.solve(A_dense, b)

    print("MNA Matrix:")
    print(A_dense)

    print("\nRHS Vector:")
    print(b)

    # Node indices follow the order in which non-ground nodes were registered.
    print("\nSolution:")
    print(f"Vin  = {x[0]:.6f} V")
    print(f"Vout = {x[1]:.6f} V")
    print(f"I_Vs = {x[2]:.6f} A")

    # Validate the computed output against the analytical divider equation.
    expected_vout = 10.0 * 2000.0 / (1000.0 + 2000.0)
    error = abs(x[1] - expected_vout)

    print("\nAnalytical check:")
    print(f"Expected Vout = {expected_vout:.6f} V")
    print(f"Error         = {error:.6e} V")


if __name__ == "__main__":
    main()
