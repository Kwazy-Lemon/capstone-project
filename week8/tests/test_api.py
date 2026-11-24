"""Exercise the public week8 API. Run: python week8/tests/test_api.py"""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit_simulator"))
from api import Capacitor, Circuit, CurrentSource, Resistor, VoltageSource


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


if __name__ == "__main__":
    unittest.main()
