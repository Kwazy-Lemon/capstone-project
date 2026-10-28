"""Zhuoyao Luo"""

import numpy as np
from scipy.sparse.linalg import spsolve

from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


def main():
    # define circuit parameters
    source_voltage = 10.0
    r1 = 1000.0
    r2 = 2000.0

    # create the circuit
    circuit = MNAAssembler()

    circuit.add_element(
        VoltageSource("Vs", "in", "gnd", source_voltage)
    )
    circuit.add_element(
        Resistor("R1", "in", "out", r1)
    )
    circuit.add_element(
        Resistor("R2", "out", "gnd", r2)
    )

    # assemble the dc equation
    g_matrix, b = circuit.assemble()

    # solve for node voltages and source current
    solution = spsolve(g_matrix, b)

    if not np.all(np.isfinite(solution)):
        raise RuntimeError("solution contains non-finite values")

    # find each unknown using its assigned index
    in_index = circuit.node_map.get_index("in")
    out_index = circuit.node_map.get_index("out")
    source_index = circuit.node_map.get_branch_index("Vs")

    vin = solution[in_index]
    vout = solution[out_index]
    source_current = solution[source_index]

    # calculate an independent analytical reference
    expected_vout = source_voltage * r2 / (r1 + r2)
    error = abs(vout - expected_vout)

    print(f"input voltage: {vin:.6f} V")
    print(f"output voltage: {vout:.6f} V")
    print(f"source current: {source_current:.6f} A")
    print(f"expected output: {expected_vout:.6f} V")
    print(f"absolute error: {error:.6e} V")

    # stop if the numerical result fails the comparison
    if not np.isclose(
        vout, expected_vout, rtol=1e-10, atol=1e-12
    ):
        raise RuntimeError("voltage divider check failed")

    print("voltage divider check passed")


if __name__ == "__main__":
    main()
