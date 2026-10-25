"""Validate new components with independently derived circuit equations.

Run from the repository root:
    python -m unittest discover -s week4/tests -v
"""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit_simulator"))
from capacitor import Capacitor
from current_source import CurrentSource
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


def node_voltage(circuit, state, node):
    index = circuit.node_map.get_index(node)
    return 0.0 if index is None else state[index]


class CapacitorTests(unittest.TestCase):
    def setUp(self):
        self.R, self.C, self.V = 1000.0, 1e-6, 5.0
        self.circuit = assemble(VoltageSource("Vin", "in", "GND", self.V),
                                Resistor("R1", "in", "out", self.R),
                                Capacitor("C1", "out", "GND", self.C))
        self.G, self.M, self.b = self.circuit.assemble_system()
        self.out = self.circuit.node_map.get_index("out")

    def initial_state(self, capacitor_voltage):
        """Satisfy the source and resistor KCL at the chosen initial voltage."""
        state = np.zeros(self.b.size)
        state[self.circuit.node_map.get_index("in")] = self.V
        state[self.out] = capacitor_voltage
        state[self.circuit.node_map.get_branch_index("Vin")] = -(self.V - capacitor_voltage) / self.R
        return state

    def test_dc_open_circuit_preserves_resistor_divider(self):
        self.circuit.add_element(Resistor("R2", "out", "GND", 2000))
        G, _, b = self.circuit.assemble_system()
        state = dc_solve(G, b)
        self.assertAlmostEqual(state[self.out], self.V * 2000 / 3000)
        # Two node voltages and one voltage-source current; C adds no unknown.
        self.assertEqual(state.size, 3)
        np.testing.assert_allclose(G @ state, b, atol=1e-13)

    def test_ac_low_pass_and_complex_source_phasor(self):
        frequencies = np.logspace(0, 6, 31)
        b_ac = self.b.astype(complex) * (1 + 2j) / self.V
        states = ac_sweep(self.G, self.M, b_ac, frequencies)
        expected = (1 + 2j) / (1 + 1j * 2 * np.pi * frequencies * self.R * self.C)
        np.testing.assert_allclose(states[:, self.out], expected, rtol=1e-12, atol=1e-14)

    def test_transient_nonzero_initial_voltage_and_discrete_recurrence(self):
        tau, dt = self.R * self.C, self.R * self.C / 20
        for initial in (0.0, 1.0, -0.5):
            with self.subTest(initial_voltage=initial):
                x0 = self.initial_state(initial)
                time, states = backward_euler(self.G, self.M, lambda t: self.b,
                                               x0, dt, 5 * tau)
                q = 1 / (1 + dt / tau)
                expected = self.V + (initial - self.V) * q ** np.arange(time.size)
                np.testing.assert_allclose(states[:, self.out], expected, atol=1e-12)
                np.testing.assert_array_equal(states[0], x0)
                # At every computed step, capacitor current equals resistor
                # current; this also checks the derivative stamp's sign.
                ic = self.C * np.diff(states[:, self.out]) / dt
                ir = (self.V - states[1:, self.out]) / self.R
                np.testing.assert_allclose(ic, ir, atol=1e-13)

    def test_transient_error_decreases_with_time_step(self):
        tau = self.R * self.C
        errors = []
        for dt in (tau / 10, tau / 20, tau / 40):
            time, states = backward_euler(self.G, self.M, lambda t: self.b,
                                           self.initial_state(0), dt, 5 * tau)
            exact = self.V * (1 - np.exp(-time / tau))
            errors.append(np.max(np.abs(states[:, self.out] - exact)))
        self.assertLess(errors[1], errors[0])
        self.assertLess(errors[2], errors[1])
        self.assertGreater(errors[1] / errors[2], 1.8)

    def test_floating_capacitor_ac_and_reversed_terminals(self):
        # A and B each have a resistor to ground; C connects A to B. Looking
        # into A gives R1 in parallel with the series combination Zc + R2.
        frequencies = np.logspace(0, 5, 21)
        r1, r2, current = 1000, 2000, 0.001
        zc = 1 / (1j * 2 * np.pi * frequencies * self.C)
        expected_a = current / (1 / r1 + 1 / (zc + r2))
        expected_b = expected_a * r2 / (zc + r2)
        for terminals in (("a", "b"), ("b", "a")):
            with self.subTest(terminals=terminals):
                circuit = assemble(Capacitor("C1", *terminals, self.C),
                                   CurrentSource("I1", "GND", "a", current),
                                   Resistor("R1", "a", "GND", r1),
                                   Resistor("R2", "b", "GND", r2))
                G, M, b = circuit.assemble_system()
                states = ac_sweep(G, M, b, frequencies)
                np.testing.assert_allclose(states[:, circuit.node_map.get_index("a")], expected_a, rtol=1e-12)
                np.testing.assert_allclose(states[:, circuit.node_map.get_index("b")], expected_b, rtol=1e-12)

    def test_floating_capacitor_transient(self):
        r1, r2, current = 1000, 2000, 0.001
        circuit = assemble(Capacitor("C1", "a", "b", self.C),
                           CurrentSource("I1", "GND", "a", current),
                           Resistor("R1", "a", "GND", r1),
                           Resistor("R2", "b", "GND", r2))
        G, M, b = circuit.assemble_system()
        a, node_b = circuit.node_map.get_index("a"), circuit.node_map.get_index("b")
        # With zero initial voltage across C, both terminals have the same
        # voltage. Their common voltage must still satisfy total resistor KCL.
        x0 = np.full(b.size, current / (1 / r1 + 1 / r2))
        tau = (r1 + r2) * self.C
        dt = tau / 30
        time, states = backward_euler(G, M, lambda t: b, x0, dt, 5 * tau)
        difference = current * r1 * (1 - (1 + dt / tau) ** (-np.arange(time.size)))
        expected_a = r1 * (current * r2 + difference) / (r1 + r2)
        expected_b = r2 * (current * r1 - difference) / (r1 + r2)
        np.testing.assert_allclose(states[:, a], expected_a, atol=1e-12)
        np.testing.assert_allclose(states[:, node_b], expected_b, atol=1e-12)

    def test_ground_at_either_capacitor_terminal(self):
        for terminals in (("out", "0"), ("gNd", "out")):
            with self.subTest(terminals=terminals):
                circuit = assemble(VoltageSource("Vin", "in", "GND", 1),
                                   Resistor("R1", "in", "out", self.R),
                                   Capacitor("C1", *terminals, self.C))
                G, M, b = circuit.assemble_system()
                state = ac_sweep(G, M, b, [1000])[0]
                expected = 1 / (1 + 1j * 2 * np.pi * 1000 * self.R * self.C)
                np.testing.assert_allclose(node_voltage(circuit, state, "out"), expected)

    def test_parallel_capacitors_add(self):
        self.circuit.add_element(Capacitor("C2", "out", "GND", 2 * self.C))
        G, M, b = self.circuit.assemble_system()
        frequencies = np.array([10, 1000, 100000])
        states = ac_sweep(G, M, b, frequencies)
        expected = self.V / (1 + 1j * 2 * np.pi * frequencies * self.R * (3 * self.C))
        np.testing.assert_allclose(states[:, self.out], expected, rtol=1e-12)

    def test_same_node_capacitor_has_no_effect(self):
        circuit = assemble(CurrentSource("I1", "GND", "out", 0.001),
                           Resistor("R1", "out", "GND", 1000),
                           Capacitor("C1", "out", "out", self.C))
        G, M, b = circuit.assemble_system()
        states = ac_sweep(G, M, b, [10, 10000])
        np.testing.assert_allclose(states[:, 0], 1.0)

    def test_invalid_capacitance(self):
        for value in (0, -1, float("nan"), float("inf"), -float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                Capacitor("C1", "a", "b", value)


class CurrentSourceTests(unittest.TestCase):
    def test_resistor_voltage_for_current_signs_and_directions(self):
        for positive, negative, direction in (("GND", "out", 1), ("out", "0", -1)):
            for current in (0.001, 0.0, -0.001):
                with self.subTest(terminals=(positive, negative), current=current):
                    circuit = assemble(CurrentSource("I1", positive, negative, current),
                                       Resistor("R1", "out", "GND", 1000))
                    G, M, b = circuit.assemble_system()
                    state = dc_solve(G, b)
                    self.assertAlmostEqual(state[0], direction * current * 1000)
                    self.assertEqual(state.size, 1)  # No source-current unknown.
                    self.assertEqual(len(circuit.branch_elements), 0)
                    self.assertEqual(M.nnz, 0)

    def test_floating_current_source_between_two_nodes(self):
        circuit = assemble(CurrentSource("I1", "a", "b", 0.001),
                           Resistor("R1", "a", "GND", 1000),
                           Resistor("R2", "b", "GND", 2000))
        G, _, b = circuit.assemble_system()
        state = dc_solve(G, b)
        self.assertAlmostEqual(node_voltage(circuit, state, "a"), -1)
        self.assertAlmostEqual(node_voltage(circuit, state, "b"), 2)
        np.testing.assert_allclose(G @ state, b, atol=1e-14)

    def test_multiple_current_sources_add_on_rhs(self):
        circuit = assemble(CurrentSource("I1", "GND", "out", 0.001),
                           CurrentSource("I2", "out", "GND", 0.00025),
                           CurrentSource("I3", "gNd", "out", 0.0005),
                           Resistor("R1", "out", "GND", 1000))
        G, _, b = circuit.assemble_system()
        state = dc_solve(G, b)
        self.assertAlmostEqual(state[0], (0.001 - 0.00025 + 0.0005) * 1000)

    def test_complex_ac_current_phasor(self):
        circuit = assemble(CurrentSource("I1", "GND", "out", 0.001),
                           Resistor("R1", "out", "GND", 1000))
        G, M, b = circuit.assemble_system()
        # The component stores a real nominal current. AC phase is specified
        # independently through the solver's complex excitation vector.
        b_ac = b.astype(complex) * (1 + 2j)
        states = ac_sweep(G, M, b_ac, [10, 1000, 100000])
        np.testing.assert_allclose(states[:, 0], 1 + 2j)

    def test_same_node_current_source_cancels(self):
        circuit = assemble(CurrentSource("I1", "out", "out", 100),
                           Resistor("R1", "out", "GND", 1000))
        G, _, b = circuit.assemble_system()
        np.testing.assert_allclose(dc_solve(G, b), 0.0)

    def test_invalid_source_current(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                CurrentSource("I1", "a", "b", value)


class CombinedCircuitTests(unittest.TestCase):
    def test_constant_current_charging_with_initial_voltage(self):
        for current in (0.001, 0.0, -0.001):
            with self.subTest(current=current):
                circuit = assemble(CurrentSource("I1", "GND", "out", current),
                                   Capacitor("C1", "out", "GND", 1e-6))
                G, M, b = circuit.assemble_system()
                time, states = backward_euler(G, M, lambda t: b, [0.25], 10e-6, 1e-3)
                expected = 0.25 + current * time / 1e-6
                np.testing.assert_allclose(states[:, 0], expected, atol=1e-12)

    def test_time_varying_current_uses_solver_callback(self):
        circuit = assemble(CurrentSource("I1", "GND", "out", 0.001),
                           Capacitor("C1", "out", "GND", 1e-6))
        G, M, b = circuit.assemble_system()
        dt, tstop = 10e-6, 1e-3
        # I(t) ramps from zero to 1 mA. The source_at callback changes only
        # the forcing, keeping component stamps and the matrix fixed.
        time, states = backward_euler(G, M, lambda t: b * t / tstop,
                                       [0.25], dt, tstop)
        n = np.arange(time.size)
        slope = 0.001 / tstop
        expected = 0.25 + slope * dt**2 / 1e-6 * n * (n + 1) / 2
        np.testing.assert_allclose(states[:, 0], expected, atol=1e-12)

    def test_current_driven_parallel_rc(self):
        current, resistance, capacitance = 0.001, 1000.0, 1e-6
        circuit = assemble(CurrentSource("I1", "GND", "out", current),
                           Resistor("R1", "out", "GND", resistance),
                           Capacitor("C1", "out", "GND", capacitance))
        G, M, b = circuit.assemble_system()
        self.assertAlmostEqual(dc_solve(G, b)[0], current * resistance)
        frequencies = np.array([10, 1000, 100000])
        states = ac_sweep(G, M, b, frequencies)
        expected = current / (1 / resistance + 1j * 2 * np.pi * frequencies * capacitance)
        np.testing.assert_allclose(states[:, 0], expected, rtol=1e-12)
        tau = resistance * capacitance
        dt = tau / 40
        time, states = backward_euler(G, M, lambda t: b, [0.2], dt, 5 * tau)
        expected = current * resistance + (0.2 - current * resistance) * (1 + dt / tau) ** (-np.arange(time.size))
        np.testing.assert_allclose(states[:, 0], expected, atol=1e-12)

    def test_parallel_rlc_with_existing_inductor(self):
        current, resistance, capacitance, inductance = 0.001, 1000, 1e-6, 0.01
        circuit = assemble(CurrentSource("I1", "GND", "out", current),
                           Resistor("R1", "out", "GND", resistance),
                           Capacitor("C1", "out", "GND", capacitance),
                           Inductor("L1", "out", "GND", inductance))
        G, M, b = circuit.assemble_system()
        dc_state = dc_solve(G, b)
        self.assertAlmostEqual(node_voltage(circuit, dc_state, "out"), 0)
        self.assertAlmostEqual(dc_state[circuit.node_map.get_branch_index("L1")], current)
        resonance = 1 / (2 * np.pi * np.sqrt(inductance * capacitance))
        frequencies = np.array([10, resonance, 100000])
        omega = 2 * np.pi * frequencies
        expected = current / (1 / resistance + 1j * omega * capacitance + 1 / (1j * omega * inductance))
        states = ac_sweep(G, M, b, frequencies)
        np.testing.assert_allclose(states[:, circuit.node_map.get_index("out")], expected, rtol=1e-12)


if __name__ == "__main__":
    unittest.main()
