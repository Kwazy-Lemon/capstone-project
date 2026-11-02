"""Jianhao Wu -- worked example of a complete, consistent RC initial state.

Place this file and initial_state.py in week5/circuit_simulator alongside
the existing components. Run: python initial_state_example.py

Circuit: 10 V source -> 1 kohm resistor -> output node -> 1 uF capacitor -> ground.
The capacitor starts at 2 V, rather than the 10 V DC equilibrium voltage.
"""

import numpy as np

from capacitor import Capacitor
from initial_state import build_initial_state, check_initial_state
from linear_solvers import backward_euler
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


def main():
    source_voltage = 10.0
    resistance = 1000.0
    capacitance = 1e-6
    initial_voltage = 2.0

    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "in", "gnd", source_voltage))
    circuit.add_element(Resistor("R1", "in", "out", resistance))
    circuit.add_element(Capacitor("C1", "out", "gnd", capacitance))

    capacitor_voltages = {"C1": initial_voltage}
    x0 = build_initial_state(circuit, capacitor_voltages)
    report = check_initial_state(circuit, x0, capacitor_voltages)

    node_in = circuit.node_map.get_index("in")
    node_out = circuit.node_map.get_index("out")
    source_index = circuit.node_map.get_branch_index("Vs")

    print("Complete initial MNA state (ground is omitted):")
    for node, index in sorted(circuit.node_map.nodes.items(), key=lambda item: item[1]):
        print(f"  x0[{index}] = V({node}) = {x0[index]:.6f} V")
    for branch, index in sorted(circuit.node_map.branches.items(), key=lambda item: item[1]):
        print(f"  x0[{index}] = I({branch}) = {x0[index]:.6f} A")
    print("  I(Vs) is positive from the source's positive to negative terminal.")

    expected_current = (source_voltage - initial_voltage) / resistance
    expected = x0.copy()
    expected[node_in] = source_voltage
    expected[node_out] = initial_voltage
    expected[source_index] = -expected_current
    if not np.allclose(x0, expected, rtol=1e-10, atol=1e-12):
        raise RuntimeError("initial state does not match the independent RC calculation")
    if not report["consistent"]:
        raise RuntimeError("valid initial state was rejected")

    print(f"\nResistor current at t=0: {expected_current:.6f} A")
    print(f"Capacitor voltage rate at t=0: {report['xdot'][node_out]:.6f} V/s")
    print("Initial-state constraint check: PASS")

    # A supplied capacitor voltage alone is not a complete MNA state.
    # Setting the source current to zero would violate KCL at the input node.
    wrong_current = x0.copy()
    wrong_current[source_index] = 0.0
    bad_report = check_initial_state(circuit, wrong_current, capacitor_voltages)
    if bad_report["consistent"]:
        raise RuntimeError("incorrect source current was not detected")
    print("\nIncorrect source-current example: correctly rejected")
    for issue in bad_report["issues"]:
        print(f"  {issue}")

    # Check the stated capacitor voltage separately from the MNA equations.
    wrong_voltage = x0.copy()
    wrong_voltage[node_out] = 3.0
    wrong_voltage[source_index] = -(source_voltage - 3.0) / resistance
    voltage_report = check_initial_state(circuit, wrong_voltage, capacitor_voltages)
    if voltage_report["consistent"]:
        raise RuntimeError("incorrect capacitor initial voltage was not detected")
    print("\nIncorrect capacitor-voltage example: correctly rejected")
    for issue in voltage_report["issues"]:
        print(f"  {issue}")

    # Reuse the existing solver, rather than implementing transient analysis again.
    # This confirms the returned vector can be passed directly to backward_euler.
    G, M, b = circuit.assemble_system()
    tau = resistance * capacitance
    dt = tau / 100
    t, X = backward_euler(G, M, lambda time: b, x0, dt, 5 * tau)
    analytical = source_voltage + (initial_voltage - source_voltage) * np.exp(-t / tau)
    error = X[:, node_out] - analytical
    if not np.array_equal(X[0], x0):
        raise RuntimeError("transient solver did not preserve the supplied t=0 state")
    if not np.allclose(X[:, node_out], analytical, rtol=0.0, atol=0.02):
        raise RuntimeError("transient result disagrees with the RC analytical solution")

    print("\nExisting backward Euler solver: PASS")
    print(f"  dt = {dt:.6e} s, tstop = {t[-1]:.6e} s")
    print(f"  Maximum sampled voltage error: {np.max(np.abs(error)):.6e} V")
    print(f"  RMS sampled voltage error: {np.sqrt(np.mean(error ** 2)):.6e} V")
    print("\nAll initial-state example checks passed.")


if __name__ == "__main__":
    main()
