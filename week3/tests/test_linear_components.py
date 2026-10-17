"""Circuit-level checks against analytical physics, not copied stamp entries.

Run from the repository root:
    python -m unittest discover -s week3/tests -v
"""

from pathlib import Path
import sys
import unittest

import numpy as np

# This repository uses standalone weekly scripts rather than an installed
# package. Resolve imports relative to this file, regardless of the shell cwd.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit_simulator"))
from controlled_sources import VCVS, VCCS, CCCS, CCVS
from inductor import Inductor
from linear_solvers import dc_solve, ac_sweep, backward_euler
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


def assemble(*elements):
    circuit = MNAAssembler()
    for element in elements:
        circuit.add_element(element)
    return circuit


def voltage(circuit, state, node):
    index = circuit.node_map.get_index(node)
    return 0.0 if index is None else state[index]


def solve_dc(circuit):
    G, _, b = circuit.assemble_system()
    x = dc_solve(G, b)
    np.testing.assert_allclose(G @ x - b, 0.0, atol=1e-12)
    return x


class RegressionTests(unittest.TestCase):
    def test_original_divider_and_dc_interface(self):
        circuit = assemble(Resistor("R1", "in", "out", 1000),
                           Resistor("R2", "out", "GND", 2000),
                           VoltageSource("Vin", "in", "GND", 10))
        G, b = circuit.assemble()  # Original two-return-value API still works.
        x = dc_solve(G, b)
        self.assertAlmostEqual(voltage(circuit, x, "out"), 10 * 2000 / 3000)

    def test_reassembly_after_adding_a_node(self):
        circuit = assemble(VoltageSource("Vin", "in", "GND", 1),
                           Inductor("L1", "in", "out", 0.01),
                           Resistor("R1", "out", "GND", 100))
        solve_dc(circuit)
        # A new node shifts all branch indices; every component must use the
        # new assignment, including a current-controlled source's reference.
        circuit.add_element(CCCS("F1", "load", "GND", "L1", 2))
        circuit.add_element(Resistor("R2", "load", "GND", 1000))
        x = solve_dc(circuit)
        self.assertAlmostEqual(voltage(circuit, x, "load"), -20.0)


class InductorTests(unittest.TestCase):
    def setUp(self):
        self.R, self.L, self.V = 100.0, 0.01, 5.0
        self.circuit = assemble(VoltageSource("Vin", "in", "GND", self.V),
                                Resistor("R1", "in", "out", self.R),
                                Inductor("L1", "out", "GND", self.L))
        self.G, self.M, self.b = self.circuit.assemble_system()
        self.k = self.circuit.node_map.get_branch_index("L1")

    def test_dc_short_circuit(self):
        x = solve_dc(self.circuit)
        self.assertAlmostEqual(voltage(self.circuit, x, "out"), 0.0)
        self.assertAlmostEqual(x[self.k], self.V / self.R)

    def test_ac_impedance(self):
        frequencies = np.logspace(0, 6, 31)
        states = ac_sweep(self.G, self.M, self.b, frequencies)
        z = 1j * 2 * np.pi * frequencies * self.L
        np.testing.assert_allclose(states[:, self.k], self.V / (self.R + z), rtol=1e-12)

    def test_transient_initial_current_and_recurrence(self):
        tau = self.L / self.R
        for initial_i in (0.0, 0.025, -0.01):
            with self.subTest(initial_current=initial_i):
                x0 = np.zeros(self.b.size)
                x0[self.circuit.node_map.get_index("in")] = self.V
                x0[self.circuit.node_map.get_index("out")] = self.V - self.R * initial_i
                x0[self.circuit.node_map.get_branch_index("Vin")] = -initial_i
                x0[self.k] = initial_i
                dt = tau / 20
                time, states = backward_euler(self.G, self.M, lambda t: self.b,
                                               x0, dt, 5 * tau)
                # General RL solution with nonzero initial current; q is the
                # backward Euler decay factor, derived from L*di/dt + R*i = V.
                final_i = self.V / self.R
                q = 1 / (1 + dt / tau)
                expected = final_i + (initial_i - final_i) * q ** np.arange(time.size)
                np.testing.assert_allclose(states[:, self.k], expected, atol=1e-13)
                np.testing.assert_array_equal(states[0], x0)

    def test_transient_error_decreases_with_time_step(self):
        x0 = np.zeros(self.b.size)
        x0[self.circuit.node_map.get_index("in")] = self.V
        x0[self.circuit.node_map.get_index("out")] = self.V
        tau = self.L / self.R
        errors = []
        for dt in (tau / 10, tau / 20, tau / 40):
            time, states = backward_euler(self.G, self.M, lambda t: self.b,
                                           x0, dt, 5 * tau)
            exact = self.V / self.R * (1 - np.exp(-time / tau))
            errors.append(np.max(np.abs(states[:, self.k] - exact)))
        self.assertLess(errors[1], errors[0])
        self.assertLess(errors[2], errors[1])
        self.assertGreater(errors[1] / errors[2], 1.8)

    def test_floating_inductor_and_reversed_direction(self):
        # Both inductor terminals are non-ground; reversing them reverses the
        # reported current, while the physical voltage across the load is unchanged.
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                terminals = ("out", "in") if reverse else ("in", "out")
                circuit = assemble(Inductor("L1", *terminals, self.L),
                                   Resistor("R1", "out", "GND", self.R),
                                   VoltageSource("Vin", "in", "GND", self.V))
                G, M, b = circuit.assemble_system()
                states = ac_sweep(G, M, b, [1000.0])
                current = self.V / (self.R + 1j * 2 * np.pi * 1000 * self.L)
                k = circuit.node_map.get_branch_index("L1")
                np.testing.assert_allclose(states[0, k], -current if reverse else current)
                np.testing.assert_allclose(voltage(circuit, states[0], "out"), self.R * current)

    def test_ground_at_positive_inductor_terminal(self):
        circuit = assemble(VoltageSource("Vin", "in", "0", self.V),
                           Resistor("R1", "in", "out", self.R),
                           Inductor("L1", "gNd", "out", self.L))
        x = solve_dc(circuit)
        self.assertAlmostEqual(x[circuit.node_map.get_branch_index("L1")], -self.V / self.R)

    def test_invalid_inductance(self):
        for value in (0, -1, float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                Inductor("L1", "a", "b", value)


class ControlledSourceTests(unittest.TestCase):
    def test_voltage_controlled_sources_with_floating_terminals_and_gains(self):
        # Vcontrol = 3 - 1 = 2 V; output negative terminal is held at 1.5 V.
        # Both control and output pairs are floating with respect to ground.
        for source_class in (VCVS, VCCS):
            for value in (0.0, 2.0, -2.0):
                with self.subTest(source=source_class.__name__, gain=value):
                    coefficient = value if source_class is VCVS else value * 0.001
                    source = source_class("D1", "out", "ref", "cp", "cm", coefficient)
                    circuit = assemble(source,
                                       VoltageSource("Vp", "cp", "GND", 3),
                                       VoltageSource("Vm", "cm", "GND", 1),
                                       VoltageSource("Vref", "ref", "GND", 1.5),
                                       Resistor("Rload", "out", "ref", 2000))
                    x = solve_dc(circuit)
                    expected_drop = coefficient * 2
                    if source_class is VCCS:
                        expected_drop *= -2000  # V/R + Iout = 0.
                    self.assertAlmostEqual(voltage(circuit, x, "out"), 1.5 + expected_drop)
                    # Ideal voltage-controlled inputs draw no current.
                    self.assertAlmostEqual(x[circuit.node_map.get_branch_index("Vp")], 0)
                    self.assertAlmostEqual(x[circuit.node_map.get_branch_index("Vm")], 0)

    def test_current_controlled_sources_with_probe_and_forward_reference(self):
        for source_class in (CCCS, CCVS):
            for value in (0.0, 2.0, -2.0):
                with self.subTest(source=source_class.__name__, gain=value):
                    coefficient = value if source_class is CCCS else value * 1000
                    source = source_class("D1", "out", "ref", "Vsense", coefficient)
                    circuit = assemble(source,
                                       Resistor("Rcontrol", "in", "sense", 1000),
                                       VoltageSource("Vsense", "sense", "GND", 0),
                                       VoltageSource("Vin", "in", "GND", 2),
                                       VoltageSource("Vref", "ref", "GND", 1.5),
                                       Resistor("Rload", "out", "ref", 2000))
                    x = solve_dc(circuit)
                    ic = x[circuit.node_map.get_branch_index("Vsense")]
                    self.assertAlmostEqual(ic, 0.002)
                    expected_drop = coefficient * 0.002
                    if source_class is CCCS:
                        expected_drop *= -2000
                    self.assertAlmostEqual(voltage(circuit, x, "out"), 1.5 + expected_drop)

    def test_grounded_and_reversed_voltage_control(self):
        for source_class, gain, expected in ((VCVS, 3.0, -6.0), (VCCS, 0.001, 4.0)):
            with self.subTest(source=source_class.__name__):
                circuit = assemble(source_class("D1", "out", "0", "gnd", "in", gain),
                                   VoltageSource("Vin", "in", "GND", 2),
                                   Resistor("Rload", "out", "GND", 2000))
                self.assertAlmostEqual(voltage(circuit, solve_dc(circuit), "out"), expected)

    def test_grounded_positive_output_terminal(self):
        sources = [
            (VCVS("E1", "GND", "out", "in", "GND", 3), -6.0),
            (VCCS("G1", "GND", "out", "in", "GND", 0.001), 4.0),
            (CCCS("F1", "GND", "out", "Vsense", 2), 8.0),
            (CCVS("H1", "GND", "out", "Vsense", 1000), -2.0),
        ]
        for source, expected in sources:
            with self.subTest(source=type(source).__name__):
                circuit = assemble(source, VoltageSource("Vin", "in", "GND", 2),
                                   Resistor("Rc", "in", "sense", 1000),
                                   VoltageSource("Vsense", "sense", "GND", 0),
                                   Resistor("Rload", "out", "GND", 2000))
                self.assertAlmostEqual(voltage(circuit, solve_dc(circuit), "out"), expected)

    def test_reversed_control_branch_current(self):
        for source_class, gain, expected in ((CCCS, 2, 8), (CCVS, 1000, -2)):
            with self.subTest(source=source_class.__name__):
                circuit = assemble(source_class("D1", "out", "GND", "Vsense", gain),
                                   VoltageSource("Vin", "in", "GND", 2),
                                   Resistor("Rc", "in", "sense", 1000),
                                   VoltageSource("Vsense", "GND", "sense", 0),
                                   Resistor("Rload", "out", "GND", 2000))
                self.assertAlmostEqual(voltage(circuit, solve_dc(circuit), "out"), expected)

    def test_current_source_can_use_inductor_control_in_ac(self):
        circuit = assemble(CCCS("F1", "load", "GND", "L1", 2),
                           VoltageSource("Vin", "in", "GND", 1),
                           Resistor("R1", "in", "out", 100),
                           Inductor("L1", "out", "GND", 0.01),
                           Resistor("Rload", "load", "GND", 1000))
        G, M, b = circuit.assemble_system()
        frequencies = np.array([10.0, 1000.0, 100000.0])
        states = ac_sweep(G, M, b, frequencies)
        expected = -2000 / (100 + 1j * 2 * np.pi * frequencies * 0.01)
        np.testing.assert_allclose(states[:, circuit.node_map.get_index("load")], expected, rtol=1e-12)

    def test_shared_output_and_control_nodes(self):
        # Vout = 2*(Vin - Vout), hence Vout = 2/3 V. Repeated node labels
        # cause several stamp entries to share a position and must be summed.
        circuit = assemble(VCVS("E1", "out", "GND", "in", "out", 2),
                           VoltageSource("Vin", "in", "GND", 1),
                           Resistor("Rload", "out", "GND", 1000))
        self.assertAlmostEqual(voltage(circuit, solve_dc(circuit), "out"), 2 / 3)
        # VCCS sensing its own voltage acts like a 1 kohm conductance.
        circuit = assemble(VCCS("G1", "out", "GND", "out", "GND", 0.001),
                           VoltageSource("Vin", "in", "GND", 1),
                           Resistor("R1", "in", "out", 1000))
        self.assertAlmostEqual(voltage(circuit, solve_dc(circuit), "out"), 0.5)

    def test_unknown_control_branch_is_rejected(self):
        for source_class in (CCCS, CCVS):
            with self.subTest(source=source_class.__name__):
                circuit = assemble(source_class("D1", "out", "GND", "missing", 1))
                with self.assertRaisesRegex(ValueError, "No branch-current unknown"):
                    circuit.assemble_system()

    def test_resistor_name_is_not_a_control_current_unknown(self):
        circuit = assemble(CCCS("F1", "out", "GND", "R1", 2),
                           Resistor("R1", "in", "GND", 1000))
        with self.assertRaisesRegex(ValueError, "No branch-current unknown"):
            circuit.assemble_system()

    def test_invalid_controlled_source_coefficients(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            for source_class in (VCVS, VCCS, CCCS, CCVS):
                with self.subTest(source=source_class.__name__, value=value):
                    with self.assertRaises(ValueError):
                        if source_class in (VCVS, VCCS):
                            source_class("D1", "out", "GND", "cp", "cm", value)
                        else:
                            source_class("D1", "out", "GND", "Vsense", value)

    def test_branch_devices_reject_stamping_before_index_assignment(self):
        for element in (Inductor("L1", "out", "GND", 0.01),
                        VCVS("E1", "out", "GND", "in", "GND", 2),
                        CCVS("H1", "out", "GND", "Vsense", 1000)):
            with self.subTest(source=type(element).__name__):
                with self.assertRaisesRegex(ValueError, "Branch index has not been assigned"):
                    element.stamp(MNAAssembler().node_map)


if __name__ == "__main__":
    unittest.main()
