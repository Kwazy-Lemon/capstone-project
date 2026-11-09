"""Zhuoyao Luo"""

from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.sparse.linalg import spsolve

# When this file is placed in week6/tests, import the simulator beside tests/.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit_simulator"))

from controlled_sources import CCVS, CCCS, VCCS, VCVS
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


class ControlledSourceTests(unittest.TestCase):
    def check_circuit(self, circuit, nodes, branches, expected_g, expected_b,
                      expected_out):
        # Assemble G x + M dx/dt = b. These DC examples have M = 0.
        g, m, b = circuit.assemble_system()
        # Matrix row/column order is all non-ground node voltages followed by
        # currents through voltage-defined branches (Vs, E1, or H1).
        self.assertEqual(circuit.node_map.nodes,
                         {name: index for index, name in enumerate(nodes)})
        self.assertEqual(circuit.node_map.branches,
                         {name: len(nodes) + index
                          for index, name in enumerate(branches)})
        # Compare against equations calculated by hand, including signs.
        # Sparse G and M are converted to arrays only for this small test.
        np.testing.assert_allclose(g.toarray(), expected_g, rtol=0, atol=1e-12)
        np.testing.assert_array_equal(m.toarray(), np.zeros_like(expected_g))
        np.testing.assert_allclose(b, expected_b, rtol=0, atol=1e-12)

        # A correct-looking stamp must also produce the expected circuit output.
        solution = spsolve(g, b)
        np.testing.assert_allclose(solution[circuit.node_map.get_index("out")],
                                   expected_out, rtol=0, atol=1e-12)

    def test_vcvs(self):
        # V(in)=1; V(out)=2 V(in). A 1 kΩ load draws 2 mA.
        # Unknowns: [V(in), V(out), I(Vs), I(E1)].
        circuit = MNAAssembler()
        circuit.add_element(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add_element(Resistor("RL", "out", "gnd", 1000.0))
        circuit.add_element(VCVS("E1", "out", "gnd", "in", "gnd", 2.0))
        # Last matrix row is V(out) - 2 V(in) = 0.
        self.check_circuit(circuit, ["in", "out"], ["Vs", "E1"],
                           [[0, 0, 1, 0], [0, .001, 0, 1],
                            [1, 0, 0, 0], [-2, 1, 0, 0]],
                           [0, 0, 1, 0], 2.0)

    def test_vccs(self):
        # 2 mS from out to ground: 2 mA leaves out, so V(out)=-2 V.
        # Unknowns: [V(in), V(out), I(Vs)]. VCCS adds no branch current.
        circuit = MNAAssembler()
        circuit.add_element(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add_element(Resistor("RL", "out", "gnd", 1000.0))
        circuit.add_element(VCCS("G1", "out", "gnd", "in", "gnd", .002))
        # Output-node KCL: 0.002 V(in) + 0.001 V(out) = 0.
        self.check_circuit(circuit, ["in", "out"], ["Vs"],
                           [[0, 0, 1], [.002, .001, 0], [1, 0, 0]],
                           [0, 0, 1], -2.0)

    def test_ccvs(self):
        # The control source supplies 1 mA, hence I(Vs)=-1 mA.
        # V(out)=1000 Ω * I(Vs)=-1 V.
        # Unknowns: [V(in), V(out), I(Vs), I(H1)].
        circuit = MNAAssembler()
        circuit.add_element(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add_element(Resistor("RC", "in", "gnd", 1000.0))
        circuit.add_element(Resistor("RL", "out", "gnd", 1000.0))
        circuit.add_element(CCVS("H1", "out", "gnd", "Vs", 1000.0))
        # Last row: V(out) - 1000 I(Vs) = 0; the gain has units of ohms.
        self.check_circuit(circuit, ["in", "out"], ["Vs", "H1"],
                           [[.001, 0, 1, 0], [0, .001, 0, 1],
                            [1, 0, 0, 0], [0, 1, -1000, 0]],
                           [0, 0, 1, 0], -1.0)

    def test_cccs(self):
        # I(Vs)=-1 mA; 2 * I(Vs) flows out to ground, so V(out)=2 V.
        # Unknowns: [V(in), V(out), I(Vs)]. CCCS adds no branch current.
        circuit = MNAAssembler()
        circuit.add_element(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add_element(Resistor("RC", "in", "gnd", 1000.0))
        circuit.add_element(Resistor("RL", "out", "gnd", 1000.0))
        circuit.add_element(CCCS("F1", "out", "gnd", "Vs", 2.0))
        # Output-node KCL: 0.001 V(out) + 2 I(Vs) = 0.
        self.check_circuit(circuit, ["in", "out"], ["Vs"],
                           [[.001, 0, 1], [0, .001, 2], [1, 0, 0]],
                           [0, 0, 1], 2.0)


if __name__ == "__main__":
    unittest.main()
