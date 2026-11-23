"""Zhuoyao Luo"""

import numpy as np

from benchmark_circuits import build_simple_rc, build_rc_ladder
from controlled_sources import VCVS
from linear_solvers import ac_sweep, dc_solve
from main import build_voltage_divider
from mna_assembler import MNAAssembler
from resistor import Resistor
from capacitor import Capacitor
from voltage_source import VoltageSource
from current_source import CurrentSource


def check_system(name, circuit, nodes, branches, expected_g, expected_m, expected_b):
    g, m, b = circuit.assemble_system()
    node_count = len(nodes)
    size = node_count + len(branches)

    # Ground is excluded; branch currents follow node voltages.
    assert circuit.node_map.num_nodes() == node_count
    assert len(circuit.branch_elements) == len(branches)
    assert circuit.node_map.nodes == {node: i for i, node in enumerate(nodes)}
    assert circuit.node_map.branches == {
        branch: node_count + i for i, branch in enumerate(branches)
    }
    assert g.shape == m.shape == (size, size)
    assert b.shape == (size,)

    # References are specified independently, without calling stamp().
    np.testing.assert_allclose(g.toarray(), expected_g, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(m.toarray(), expected_m, rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(b, expected_b, rtol=1e-12, atol=1e-15)
    print(f"{name}: PASS (nodes={node_count}, branches={len(branches)}, size={size})")
    return g, m, b


def check_components():
    # Floating terminals exercise both diagonal and off-diagonal entries.
    circuit = MNAAssembler()
    circuit.add_element(Resistor("R1", "p", "n", 1000.0))
    check_system(
        "Resistor", circuit, ["p", "n"], [],
        [[0.001, -0.001], [-0.001, 0.001]],
        np.zeros((2, 2)), [0, 0],
    )

    circuit = MNAAssembler()
    circuit.add_element(Capacitor("C1", "p", "n", 1e-6))
    check_system(
        "Capacitor", circuit, ["p", "n"], [],
        np.zeros((2, 2)),
        [[1e-6, -1e-6], [-1e-6, 1e-6]], [0, 0],
    )

    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "p", "n", 5.0))
    check_system(
        "Voltage source", circuit, ["p", "n"], ["Vs"],
        [[0, 0, 1], [0, 0, -1], [1, -1, 0]],
        np.zeros((3, 3)), [0, 0, 5],
    )


def check_benchmarks():
    conductance = 1 / 1000.0
    capacitance = 1e-6

    check_system(
        "Voltage divider", build_voltage_divider(10.0, 1000.0, 2000.0),
        ["in", "out"], ["Vs"],
        [[0.001, -0.001, 1],
         [-0.001, 0.0015, 0],
         [1, 0, 0]],
        np.zeros((3, 3)), [0, 0, 10],
    )

    # Unknown order: V(in), V(out), I(Vs).
    check_system(
        "Simple RC", build_simple_rc(), ["in", "out"], ["Vs"],
        [[conductance, -conductance, 1],
         [-conductance, conductance, 0],
         [1, 0, 0]],
        np.diag([0, capacitance, 0]), [0, 0, 1],
    )

    # Unknown order: V(in), V(n1), ..., V(n4), I(Vs).
    g = conductance
    check_system(
        "Four-stage RC", build_rc_ladder(),
        ["in", "n1", "n2", "n3", "n4"], ["Vs"],
        [[g, -g, 0, 0, 0, 1],
         [-g, 2*g, -g, 0, 0, 0],
         [0, -g, 2*g, -g, 0, 0],
         [0, 0, -g, 2*g, -g, 0],
         [0, 0, 0, -g, g, 0],
         [1, 0, 0, 0, 0, 0]],
        np.diag([0, capacitance, capacitance, capacitance, capacitance, 0]),
        [0, 0, 0, 0, 0, 1],
    )


def check_vcvs_feedback():
    # Buffer the RC node with a finite-gain VCVS and direct negative feedback.
    resistance = 1000.0
    capacitance = 1e-6
    gain = 1000.0
    conductance = 1 / resistance
    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "in", "0", 1.0))
    circuit.add_element(Resistor("R1", "in", "rc", resistance))
    circuit.add_element(Capacitor("C1", "rc", "0", capacitance))
    circuit.add_element(VCVS("E1", "out", "0", "rc", "out", gain))

    # V(out) - gain * (V(rc) - V(out)) = 0.
    g, m, b = check_system(
        "VCVS negative-feedback RC low-pass", circuit,
        ["in", "rc", "out"], ["Vs", "E1"],
        [[conductance, -conductance, 0, 1, 0],
         [-conductance, conductance, 0, 0, 0],
         [0, 0, 0, 0, 1],
         [1, 0, 0, 0, 0],
         [0, -gain, 1 + gain, 0, 0]],
        np.diag([0, capacitance, 0, 0, 0]), [0, 0, 0, 1, 0],
    )

    input_index = circuit.node_map.get_index("in")
    output_index = circuit.node_map.get_index("out")
    dc = dc_solve(g, b)
    closed_loop_gain = gain / (1 + gain)
    np.testing.assert_allclose(
        dc[output_index] / dc[input_index], closed_loop_gain, rtol=1e-12,
    )

    cutoff_hz = 1 / (2 * np.pi * resistance * capacitance)
    ac_source = np.array([0, 0, 0, 1, 0], dtype=complex)
    ac = ac_sweep(g, m, ac_source, [cutoff_hz])
    np.testing.assert_allclose(
        ac[0, output_index] / ac[0, input_index],
        closed_loop_gain / (1 + 1j), rtol=1e-12,
    )
    print("VCVS closed-loop DC gain and cutoff response: PASS")


def check_addition():
    circuit = MNAAssembler()
    circuit.add_element(Resistor("R1", "out", "gnd", 1000.0))
    circuit.add_element(Resistor("R2", "gnd", "out", 2000.0))
    circuit.add_element(Capacitor("C1", "out", "gnd", 1e-6))
    circuit.add_element(Capacitor("C2", "gnd", "out", 2e-6))
    circuit.add_element(CurrentSource("I1", "gnd", "out", 0.001))
    circuit.add_element(CurrentSource("I2", "out", "gnd", 0.0004))

    # G = 1/R1 + 1/R2; M = C1 + C2; b = I1 - I2.
    check_system(
        "Shared entries", circuit, ["out"], [],
        [[0.0015]], [[3e-6]], [0.0006],
    )


if __name__ == "__main__":
    check_components()
    check_benchmarks()
    check_vcvs_feedback()
    check_addition()
    print("All assembly checks passed.")
