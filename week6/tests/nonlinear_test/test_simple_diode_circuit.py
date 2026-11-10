"""Regression test for the simple diode circuit's DC operating point."""

from pathlib import Path
import sys
import unittest

# The IDE can load this test from the project root.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from simple_diode_circuit import main


class SimpleDiodeCircuitTest(unittest.TestCase):
    def test_dc_operating_point(self):
        # The example raises if the solver fails or the circuit checks fail.
        main()


if __name__ == "__main__":
    unittest.main()
