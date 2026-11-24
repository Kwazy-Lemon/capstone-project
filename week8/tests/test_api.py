"""Exercise the public week8 API. Run: python week8/tests/test_api.py"""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit_simulator"))
from api import (
    Capacitor, CCCS, CCVS, Circuit, CurrentSource, Inductor,
    Resistor, Result, VCCS, VCVS, VoltageSource,
)


class TestCircuitAPI(unittest.TestCase):
    def test_dc_voltage_divider(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "gnd", 10.0))
        circuit.add(Resistor("R1", "in", "out", 1000.0))
        circuit.add(Resistor("R2", "out", "gnd", 2000.0))

        result = circuit.dc()
        self.assertAlmostEqual(result.voltage("out"), 20 / 3)
        self.assertAlmostEqual(result.current("Vs"), -10 / 3000)
        self.assertEqual(result.voltage("0"), 0.0)

    def test_rc_low_pass(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "gnd", 5.0))
        circuit.add(Resistor("R1", "in", "out", 1000.0))
        circuit.add(Capacitor("C1", "out", "gnd", 1e-6))

        cutoff = 1 / (2 * np.pi * 1000 * 1e-6)
        frequencies = np.array([10.0, cutoff, 1000.0])
        result = circuit.ac(frequencies, sources={"Vs": 1.0})
        expected = 1 / (1 + 1j * frequencies / cutoff)
        np.testing.assert_allclose(result.voltage("in"), np.ones(3))
        np.testing.assert_allclose(result.voltage("out"), expected)
        np.testing.assert_allclose(result.current("Vs"), (expected - 1) / 1000)
        np.testing.assert_array_equal(result.voltage("gnd"), np.zeros(3))
        self.assertAlmostEqual(circuit.dc().voltage("out"), 5.0)

    def test_ac_current_source_direction(self):
        circuit = Circuit()
        circuit.add(Resistor("R1", "a", "gnd", 1000.0))
        circuit.add(Resistor("R2", "b", "gnd", 2000.0))
        circuit.add(CurrentSource("Is", "a", "b", 0.0))

        result = circuit.ac([100.0], sources={"Is": 1e-3j})
        np.testing.assert_allclose(result.voltage("a"), [-1j])
        np.testing.assert_allclose(result.voltage("b"), [2j])

    def test_loaded_four_stage_rc_ladder(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "gnd", 5.0))
        resistance = 1000.0
        capacitance = 1e-6
        load = 2000.0
        nodes = ["in", "n1", "n2", "n3", "out"]
        for stage in range(1, len(nodes)):
            circuit.add(Resistor(f"R{stage}", nodes[stage - 1], nodes[stage], resistance))
            circuit.add(Capacitor(f"C{stage}", nodes[stage], "gnd", capacitance))
        circuit.add(Resistor("RL", "out", "gnd", load))

        # At DC all capacitors are open and the same current flows through all resistors.
        dc = circuit.dc()
        current = 5.0 / (4 * resistance + load)
        for stage, node in enumerate(nodes):
            self.assertAlmostEqual(dc.voltage(node), 5.0 - stage * resistance * current)
        self.assertAlmostEqual(dc.current("Vs"), -current)

        frequencies = np.array([10.0, 100.0, 1000.0])
        ac = circuit.ac(frequencies, sources={"Vs": 1.0})
        # Independent AC reference: assume Vout = 1 V, work back to Vin using
        # KCL and Ohm's law, then scale all voltages to the actual 1 V input.
        y_cap = 2j * np.pi * frequencies * capacitance
        voltage = np.ones(len(frequencies), dtype=complex)
        downstream_current = voltage / load
        reference = {"out": voltage}
        for stage in range(4, 0, -1):
            downstream_current = downstream_current + y_cap * voltage
            voltage = voltage + resistance * downstream_current
            reference[nodes[stage - 1]] = voltage
        for node in nodes:
            np.testing.assert_allclose(ac.voltage(node), reference[node] / reference["in"])

    def test_vcvs_active_low_pass(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "gnd", 0.5))
        circuit.add(Resistor("R1", "in", "sense", 1000.0))
        circuit.add(Capacitor("C1", "sense", "gnd", 1e-6))
        circuit.add(VCVS("E1", "out", "gnd", "sense", "gnd", 4.0))
        circuit.add(Resistor("RL", "out", "gnd", 2000.0))

        dc = circuit.dc()
        self.assertAlmostEqual(dc.voltage("out"), 2.0)
        self.assertAlmostEqual(dc.current("E1"), -1e-3)

        cutoff = 1 / (2 * np.pi * 1000 * 1e-6)
        frequencies = cutoff * np.array([0.1, 1.0, 10.0])
        ac = circuit.ac(frequencies, sources={"Vs": 1.0})
        expected = 4 / (1 + 1j * frequencies / cutoff)
        np.testing.assert_allclose(ac.voltage("out"), expected)
        np.testing.assert_allclose(ac.current("E1"), -expected / 2000.0)

    def test_ac_source_superposition(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "gnd", 5.0))
        circuit.add(Resistor("R1", "in", "out", 1000.0))
        circuit.add(Resistor("R2", "out", "gnd", 1000.0))
        circuit.add(CurrentSource("Is", "gnd", "out", 0.0))

        # At out, KCL gives Vout = Vs/2 + (500 ohms)*Is.
        for sources, expected in [
            ({"Vs": 2.0}, 1.0),
            ({"Is": 2e-3j}, 1j),
            ({"Vs": 2.0, "Is": 2e-3j}, 1 + 1j),
        ]:
            with self.subTest(sources=sources):
                ac = circuit.ac([100.0], sources)
                np.testing.assert_allclose(ac.voltage("out"), [expected])
        self.assertAlmostEqual(circuit.dc().voltage("out"), 2.5)

    def test_invalid_ac_source(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "gnd", 1.0))
        circuit.add(Resistor("R1", "in", "gnd", 1000.0))
        for name in ["missing", "R1"]:
            with self.subTest(name=name), self.assertRaises(ValueError):
                circuit.ac([100.0], sources={name: 1.0})

    def test_inductor_dc_and_ac(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "gnd", 2.0))
        circuit.add(Resistor("R1", "in", "out", 100.0))
        circuit.add(Inductor("L1", "out", "gnd", 0.01))

        dc = circuit.dc()
        self.assertIsInstance(dc, Result)
        self.assertAlmostEqual(dc.voltage("out"), 0.0)
        self.assertAlmostEqual(dc.current("L1"), 0.02)

        # Preserve the requested frequency order, including repeated points.
        frequencies = np.array([1000.0, 10.0, 1000.0])
        z_l = 2j * np.pi * frequencies * 0.01
        current = (1 + 2j) / (100 + z_l)
        ac = circuit.ac(frequencies, {"Vs": 1 + 2j})
        self.assertIsInstance(ac, Result)
        np.testing.assert_allclose(ac.current("L1"), current)
        np.testing.assert_allclose(ac.voltage("out"), z_l * current)

    def test_vccs_differential_control(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vp", "p", "0", 3.0))
        circuit.add(VoltageSource("Vn", "n", "0", 1.0))
        circuit.add(VCCS("G1", "out", "0", "p", "n", 0.002))
        circuit.add(Resistor("RL", "out", "0", 1000.0))

        # Vout = -RL * gm * (Vp - Vn).
        self.assertAlmostEqual(circuit.dc().voltage("out"), -4.0)
        ac = circuit.ac([100], {"Vp": 1 + 1j, "Vn": 2j})
        np.testing.assert_allclose(ac.voltage("out"), [-2 + 2j])

    def test_cccs_branch_control(self):
        circuit = Circuit()
        # Refer to the controlling branch before adding it.
        circuit.add(CCCS("F1", "out", "0", "Vs", 3.0))
        circuit.add(Resistor("RL", "out", "0", 2000.0))
        circuit.add(VoltageSource("Vs", "in", "0", 2.0))
        circuit.add(Resistor("Rin", "in", "0", 1000.0))

        # I(Vs) = -Vin / Rin and Vout = -RL * 3 * I(Vs).
        dc = circuit.dc()
        self.assertAlmostEqual(dc.current("Vs"), -0.002)
        self.assertAlmostEqual(dc.voltage("out"), 12.0)
        ac = circuit.ac([100], {"Vs": 1j})
        np.testing.assert_allclose(ac.voltage("out"), [6j])

    def test_ccvs_branch_control(self):
        circuit = Circuit()
        circuit.add(CCVS("H1", "out", "0", "Vs", 500.0))
        circuit.add(Resistor("RL", "out", "0", 2000.0))
        circuit.add(VoltageSource("Vs", "in", "0", 2.0))
        circuit.add(Resistor("Rin", "in", "0", 1000.0))

        # Vout = 500 ohms * I(Vs) and I(H1) = -Vout / RL.
        dc = circuit.dc()
        self.assertAlmostEqual(dc.voltage("out"), -1.0)
        self.assertAlmostEqual(dc.current("H1"), 0.0005)
        ac = circuit.ac([100], {"Vs": 1j})
        np.testing.assert_allclose(ac.voltage("out"), [-0.5j])
        np.testing.assert_allclose(ac.current("H1"), [0.00025j])

    def test_dc_current_source_and_unexcited_ac(self):
        circuit = Circuit()
        circuit.add(CurrentSource("Is", "GND", "out", 0.003))
        circuit.add(Resistor("R1", "out", "0", 2000.0))
        self.assertAlmostEqual(circuit.dc().voltage("out"), 6.0)
        ac = circuit.ac([10, 20], {})
        np.testing.assert_array_equal(ac.voltage("out"), [0.0, 0.0])
        np.testing.assert_array_equal(ac.voltage("GND"), [0.0, 0.0])

    def test_saved_results_after_adding_components(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "0", 1.0))
        circuit.add(Resistor("R1", "in", "0", 1000.0))
        before = circuit.dc()
        # A new node moves the source-current index when the circuit is reassembled.
        circuit.add(Resistor("R2", "in", "new", 500.0))
        circuit.add(Resistor("R3", "new", "0", 500.0))
        after = circuit.dc()
        self.assertAlmostEqual(after.voltage("new"), 0.5)
        self.assertAlmostEqual(after.current("Vs"), -0.002)
        self.assertAlmostEqual(before.current("Vs"), -0.001)
        with self.assertRaises(KeyError):
            before.voltage("new")

    def test_invalid_result_names(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "0", 1.0))
        circuit.add(Resistor("R1", "in", "0", 1000.0))
        for result in [circuit.dc(), circuit.ac([100], {"Vs": 1.0})]:
            with self.assertRaises(KeyError):
                result.voltage("missing")
            with self.assertRaises(KeyError):
                result.current("R1")  # No MNA branch current for a resistor.

    def test_invalid_sweep_inputs(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "0", 1.0))
        for frequencies in [[], [0], [-1], [np.nan], [np.inf], [[100]]]:
            with self.subTest(frequencies=frequencies), self.assertRaises(ValueError):
                circuit.ac(frequencies, {"Vs": 1.0})
        for phasor in [np.nan, np.inf]:
            with self.subTest(phasor=phasor), self.assertRaises(ValueError):
                circuit.ac([100], {"Vs": phasor})

    def test_duplicate_name_preserves_circuit(self):
        circuit = Circuit()
        circuit.add(VoltageSource("Vs", "in", "0", 1.0))
        circuit.add(Resistor("R1", "in", "0", 1000.0))
        with self.assertRaises(ValueError):
            circuit.add(Resistor("R1", "in", "0", 500.0))
        self.assertAlmostEqual(circuit.dc().current("Vs"), -0.001)

    def test_missing_control_branch(self):
        for element in [CCCS("F1", "out", "0", "missing", 2.0),
                        CCVS("H1", "out", "0", "missing", 500.0)]:
            with self.subTest(source=element.name):
                circuit = Circuit()
                circuit.add(Resistor("RL", "out", "0", 1000.0))
                circuit.add(element)
                with self.assertRaises(ValueError):
                    circuit.dc()


if __name__ == "__main__":
    unittest.main()
