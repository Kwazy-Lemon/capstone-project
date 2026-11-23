"""Run the Week 2 MATLAB resistor–diode circuit with the Python DC solver."""

from math import isclose
from pathlib import Path
import sys

import numpy as np

# Use the Week 6 simulator modules when running from the tests folder.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "circuit_simulator"))

from mna_assembler import MNAAssembler
from nonlinear.diode import Diode
from nonlinear.nonlinear_solver import NonlinearSolver
from resistor import Resistor
from voltage_source import VoltageSource


def main():
    supply_voltage = 5.0
    resistance = 1000.0
    saturation_current = 1e-9
    thermal_voltage = 0.02585
    ideality_factor = 1.0
    initial_diode_voltage = 0.4
    voltage_tol = 1e-9
    current_tol = 1e-12

    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "in", "0", supply_voltage))
    circuit.add_element(Resistor("R1", "in", "out", resistance))
    circuit.add_element(
        Diode(
            "D1", "out", "0",
            Is=saturation_current, n=ideality_factor, Vt=thermal_voltage,
        )
    )

    # MNA ordering for this circuit: input voltage, diode voltage, source current.
    initial_source_current = -(supply_voltage - initial_diode_voltage) / resistance
    x0 = [supply_voltage, initial_diode_voltage, initial_source_current]
    x = NonlinearSolver(circuit).solve(
        max_iter=100, tol=voltage_tol, residual_tol=current_tol, x0=x0
    )

    v_in = x[circuit.node_map.get_index("in")]
    v_diode = x[circuit.node_map.get_index("out")]
    i_source = x[circuit.node_map.get_branch_index("Vs")]
    i_resistor = (v_in - v_diode) / resistance
    i_diode = saturation_current * np.expm1(
        v_diode / (ideality_factor * thermal_voltage)
    )

    # Check the source constraint and both node current balances.
    if not (
        isclose(v_in, supply_voltage, rel_tol=0.0, abs_tol=voltage_tol)
        and isclose(i_resistor, -i_source, rel_tol=0.0, abs_tol=current_tol)
        and isclose(i_diode, i_resistor, rel_tol=0.0, abs_tol=current_tol)
    ):
        raise RuntimeError("The operating point does not satisfy the circuit equations.")

    print("Converged: simple diode circuit")
    print(f"Source voltage: {v_in:.9f} V")
    print(f"Diode voltage:  {v_diode:.9f} V")
    print(f"Diode current:  {i_diode * 1e3:.9f} mA")
    print(f"Source current: {i_source * 1e3:.9f} mA")


if __name__ == "__main__":
    main()
