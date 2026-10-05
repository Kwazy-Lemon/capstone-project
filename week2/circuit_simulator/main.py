import numpy as np

from resistor import Resistor
from voltage_source import VoltageSource
from mna_assembler import MNAAssembler


def main():
    # Create the circuit
    assembler = MNAAssembler()

    # Circuit values
    R1 = Resistor("R1", "Vin", "out", 1000.0)
    R2 = Resistor("R2", "out", "GND", 2000.0)
    Vs = VoltageSource("Vs", "Vin", "GND", 10.0)

    # Add elements to the circuit
    assembler.add_element(R1)
    assembler.add_element(R2)
    assembler.add_element(Vs)

    # Assemble the MNA system
    A, b = assembler.assemble()

    # Convert sparse matrix to dense format for display
    A_dense = A.toarray()

    # Solve Ax = b
    x = np.linalg.solve(A_dense, b)

    # Display results
    print("MNA Matrix:")
    print(A_dense)

    print("\nRHS Vector:")
    print(b)

    print("\nSolution:")
    print(f"Vin  = {x[0]:.6f} V")
    print(f"Vout = {x[1]:.6f} V")
    print(f"I_Vs = {x[2]:.6f} A")

    # Analytical check
    expected_vout = 10.0 * 2000.0 / (1000.0 + 2000.0)

    print("\nAnalytical check:")
    print(f"Expected Vout = {expected_vout:.6f} V")
    print(f"Error         = {abs(x[1] - expected_vout):.6e} V")


if __name__ == "__main__":
    main()
