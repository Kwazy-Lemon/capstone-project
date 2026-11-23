"""Zhuoyao Luo"""

import unittest

import numpy as np
from scipy.sparse.linalg import spsolve

from controlled_sources import CCVS, CCCS, VCCS, VCVS
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


class ControlledSourceTests(unittest.TestCase):
    def check_circuit(self, circuit, nodes, branches, expected_g, expected_b,
                      expected_out):
        # Assemble G x + M dx/dt = b. These DC examples have M = 0.
        g, m, b = circuit.assemble_system()

        # Node voltages come first; branch currents come after them.
        self.assertEqual(
            circuit.node_map.nodes,
            {name: index for index, name in enumerate(nodes)}
        )
        self.assertEqual(
            circuit.node_map.branches,
            {
                name: len(nodes) + index
                for index, name in enumerate(branches)
            }
        )

        # Compare the assembled system with independently calculated values.
        np.testing.assert_allclose(
            g.toarray(), expected_g, rtol=0, atol=1e-12
        )
        np.testing.assert_array_equal(
            m.toarray(), np.zeros_like(expected_g)
        )
        np.testing.assert_allclose(
            b, expected_b, rtol=0, atol=1e-12
        )

        # Solve the circuit and check its output voltage.
        solution = spsolve(g, b)
        np.testing.assert_allclose(
            solution[circuit.node_map.get_index("out")],
            expected_out,
            rtol=0,
            atol=1e-12
        )

    def test_vcvs(self):
        # V(in) = 1 V; V(out) = 2 V(in) = 2 V.
        # Unknowns: [V(in), V(out), I(Vs), I(E1)].
        circuit = MNAAssembler()
        circuit.add_element(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add_element(Resistor("RL", "out", "gnd", 1000.0))
        circuit.add_element(
            VCVS("E1", "out", "gnd", "in", "gnd", 2.0)
        )

        # The last matrix row represents V(out) - 2 V(in) = 0.
        self.check_circuit(
            circuit,
            ["in", "out"],
            ["Vs", "E1"],
            [
                [0, 0, 1, 0],
                [0, .001, 0, 1],
                [1, 0, 0, 0],
                [-2, 1, 0, 0],
            ],
            [0, 0, 1, 0],
            2.0,
        )

    def test_vccs(self):
        # The source sends 2 mA from out to ground.
        # With a 1 kΩ load, V(out) = -2 V.
        # Unknowns: [V(in), V(out), I(Vs)].
        circuit = MNAAssembler()
        circuit.add_element(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add_element(Resistor("RL", "out", "gnd", 1000.0))
        circuit.add_element(
            VCCS("G1", "out", "gnd", "in", "gnd", .002)
        )

        # Output-node KCL: 0.002 V(in) + 0.001 V(out) = 0.
        self.check_circuit(
            circuit,
            ["in", "out"],
            ["Vs"],
            [
                [0, 0, 1],
                [.002, .001, 0],
                [1, 0, 0],
            ],
            [0, 0, 1],
            -2.0,
        )

    def test_ccvs(self):
        # The voltage source supplies 1 mA, so I(Vs) = -1 mA
        # in the MNA reference direction.
        # V(out) = 1000 Ω × I(Vs) = -1 V.
        # Unknowns: [V(in), V(out), I(Vs), I(H1)].
        circuit = MNAAssembler()
        circuit.add_element(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add_element(Resistor("RC", "in", "gnd", 1000.0))
        circuit.add_element(Resistor("RL", "out", "gnd", 1000.0))
        circuit.add_element(
            CCVS("H1", "out", "gnd", "Vs", 1000.0)
        )

        # The last row represents V(out) - 1000 I(Vs) = 0.
        self.check_circuit(
            circuit,
            ["in", "out"],
            ["Vs", "H1"],
            [
                [.001, 0, 1, 0],
                [0, .001, 0, 1],
                [1, 0, 0, 0],
                [0, 1, -1000, 0],
            ],
            [0, 0, 1, 0],
            -1.0,
        )

    def test_cccs(self):
        # I(Vs) = -1 mA. The CCCS output current is 2 I(Vs),
        # which gives V(out) = 2 V across the 1 kΩ load.
        # Unknowns: [V(in), V(out), I(Vs)].
        circuit = MNAAssembler()
        circuit.add_element(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add_element(Resistor("RC", "in", "gnd", 1000.0))
        circuit.add_element(Resistor("RL", "out", "gnd", 1000.0))
        circuit.add_element(
            CCCS("F1", "out", "gnd", "Vs", 2.0)
        )

        # Output-node KCL: 0.001 V(out) + 2 I(Vs) = 0.
        self.check_circuit(
            circuit,
            ["in", "out"],
            ["Vs"],
            [
                [.001, 0, 1],
                [0, .001, 2],
                [1, 0, 0],
            ],
            [0, 0, 1],
            2.0,
        )


if __name__ == "__main__":
    unittest.main()
