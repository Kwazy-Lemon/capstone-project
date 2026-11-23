"""Zhuoyao Luo"""

from mna_assembler import MNAAssembler
from resistor import Resistor
from capacitor import Capacitor
from voltage_source import VoltageSource


def build_simple_rc(voltage=1.0):
    """R1: in -> out, 1 kohm. C1: out -> gnd, 1 uF."""
    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "in", "gnd", voltage))
    circuit.add_element(Resistor("R1", "in", "out", 1000.0))
    circuit.add_element(Capacitor("C1", "out", "gnd", 1e-6))
    return circuit


def build_rc_ladder(voltage=1.0):
    """Four stages: in -> n1 -> n2 -> n3 -> n4 (output).

    Each series resistor is 1 kohm. Each stage node has 1 uF to ground.
    """
    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "in", "gnd", voltage))

    previous = "in"
    for stage in range(1, 5):
        node = f"n{stage}"
        circuit.add_element(Resistor(f"R{stage}", previous, node, 1000.0))
        circuit.add_element(Capacitor(f"C{stage}", node, "gnd", 1e-6))
        previous = node

    return circuit
