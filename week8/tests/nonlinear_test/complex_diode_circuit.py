"""DC-only network with series, parallel, and reverse-biased diode branches."""

import csv
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "circuit_simulator"))

from current_source import CurrentSource
from mna_assembler import MNAAssembler
from nonlinear.diode import Diode
from nonlinear.nonlinear_solver import NonlinearSolver
from resistor import Resistor
from voltage_source import VoltageSource


def main():
    saturation_current = 1e-9
    ideality_factor = 1.0
    # Match .temp 27 and Tnom=27 in the LTspice model.
    thermal_voltage = 0.0258649529170901
    voltage_tol = 1e-9
    current_tol = 1e-12

    circuit = MNAAssembler()
    voltage_sources = [
        VoltageSource("V1", "vp", "0", 5.0),
        VoltageSource("V2", "vb", "0", 2.0),
    ]
    current_sources = [
        CurrentSource("I1", "0", "c", 0.0004),
        CurrentSource("I2", "d", "0", 0.0001),
    ]
    resistors = [
        Resistor("R1", "vp", "a", 1000.0),
        Resistor("R2", "vb", "b", 2200.0),
        Resistor("R3", "a", "b", 3300.0),
        Resistor("R4", "a", "c", 1500.0),
        Resistor("R5", "b", "d", 1000.0),
        Resistor("R6", "c", "d", 4700.0),
        Resistor("R7", "c", "0", 10000.0),
        Resistor("R8", "d", "0", 6800.0),
        Resistor("R9", "b", "0", 4700.0),
        Resistor("R10", "a", "d", 10000.0),
    ]
    # D1/D2 share both terminals. Each series pair has an unbranched midpoint.
    diode_connections = [
        ("D1", "a", "0"), ("D2", "a", "0"),
        ("D3", "b", "e"), ("D4", "e", "0"),
        ("D5", "c", "f"), ("D6", "f", "0"),
        ("D7", "b", "c"), ("D8", "d", "0"),
        ("D9", "0", "a"), ("D10", "c", "b"),
    ]
    diodes = [
        Diode(name, anode, cathode, Is=saturation_current,
              n=ideality_factor, Vt=thermal_voltage)
        for name, anode, cathode in diode_connections
    ]
    for element in voltage_sources + current_sources + resistors + diodes:
        circuit.add_element(element)

    # A zero initial guess exercises the existing full-step Newton solver.
    solver = NonlinearSolver(circuit)
    with np.errstate(over="raise", invalid="raise"):
        x = solver.solve_dc(max_iter=1000, tol=voltage_tol, residual_tol=current_tol)
    G, _, b = circuit.assemble_system()
    residual = solver.get_residual(x, G, b)
    np.testing.assert_allclose(residual[:circuit.node_map.num_nodes()], 0,
                               rtol=0, atol=current_tol, err_msg="Node KCL")
    np.testing.assert_allclose(residual[circuit.node_map.num_nodes():], 0,
                               rtol=0, atol=voltage_tol, err_msg="Voltage sources")

    voltages = {"0": 0.0}
    for node, index in circuit.node_map.nodes.items():
        voltages[node] = x[index]
    results = {f"V({node})": value for node, value in voltages.items() if node != "0"}
    power = 0.0
    for element in resistors:
        drop = voltages[element.node1] - voltages[element.node2]
        current = drop / element.resistance
        results[f"I({element.name})"] = current
        power += drop * current
    diode_currents = {}
    for diode in diodes:
        drop = voltages[diode.anode] - voltages[diode.cathode]
        current = diode.Is * np.expm1(drop / (diode.n * diode.Vt))
        diode_currents[diode.name] = current
        results[f"I({diode.name})"] = current
        power += drop * current
    for source in voltage_sources:
        current = x[circuit.node_map.get_branch_index(source.name)]
        results[f"I({source.name})"] = current
        power += (voltages[source.node_plus] - voltages[source.node_minus]) * current
    for source in current_sources:
        results[f"I({source.name})"] = source.current
        power += (voltages[source.node_plus] - voltages[source.node_minus]) * source.current

    for first, second in (("D1", "D2"), ("D3", "D4"), ("D5", "D6")):
        np.testing.assert_allclose(diode_currents[first], diode_currents[second],
                                   rtol=0, atol=current_tol)
    if diode_currents["D9"] >= 0 or diode_currents["D10"] >= 0:
        raise AssertionError("Expected D9 and D10 to be reverse-biased.")
    np.testing.assert_allclose(power, 0, rtol=0, atol=1e-10, err_msg="Power balance")

    # Captured from a real LTspice 26.0.2 .op run of the matching schematic.
    reference_path = Path(__file__).with_name("complex_diode_ltspice_reference.csv")
    with reference_path.open(newline="", encoding="utf-8") as reference_file:
        reference = {
            row["quantity"]: float(row["ltspice_value"])
            for row in csv.DictReader(reference_file)
        }
    if results.keys() != reference.keys():
        raise AssertionError("Python and LTspice must report the same quantities.")

    print("Converged: complex diode circuit (DC operating point)")
    print("Current is positive from the first terminal to the second.")
    print(f"{'Quantity':<12} {'Python':>16} {'LTspice':>16} {'Abs. error':>13} {'Unit':>5}")
    for name, value in results.items():
        expected = reference[name]
        unit = "V" if name.startswith("V(") else "A"
        print(f"{name:<12} {value:16.9g} {expected:16.9g} {abs(value-expected):13.3e} {unit:>5}")
        comparison_tol = 1e-7 if unit == "V" else 2e-12
        np.testing.assert_allclose(value, expected, rtol=1e-6, atol=comparison_tol,
                                   err_msg=f"LTspice comparison for {name}")
    print(f"PASS: all {len(results)} values match the saved LTspice operating point.")
    #print(f"Maximum KCL residual: {np.max(np.abs(residual[:circuit.node_map.num_nodes()])):.3e} A")
    #print(f"Power balance residual: {power:.3e} W")
    return results


if __name__ == "__main__":
    main()
