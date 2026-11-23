"""Zhuoyao Luo"""

import numpy as np

from benchmark_circuits import build_simple_rc, build_rc_ladder
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
    check_addition()
    print("All assembly checks passed.")
