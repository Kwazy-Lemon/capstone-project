"""Public Python interface to the week8 linear DC and AC solvers."""

import numpy as np

from capacitor import Capacitor
from controlled_sources import CCCS, CCVS, VCCS, VCVS
from current_source import CurrentSource
from inductor import Inductor
from linear_solvers import ac_sweep, dc_solve
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource

__all__ = [
    "Circuit", "Result", "Resistor", "Capacitor", "Inductor",
    "VoltageSource", "CurrentSource", "VCVS", "VCCS", "CCVS", "CCCS",
]


class Result:
    """Named DC values or AC traces in the requested frequency order."""

    def __init__(self, values, node_map):
        self._values = values
        self._nodes = node_map.nodes.copy()
        self._branches = node_map.branches.copy()

    def voltage(self, node):
        """Return the voltage at a node relative to ground, in volts."""
        if node.lower() in ("0", "gnd"):
            if self._values.ndim == 1:
                return 0.0
            return np.zeros(self._values.shape[0], dtype=self._values.dtype)
        return self._values[..., self._nodes[node]]

    def current(self, branch):
        """Return an MNA branch current, e.g. the current through Vs, in amps."""
        return self._values[..., self._branches[branch]]


class Circuit:
    def __init__(self):
        self._mna = MNAAssembler()

    def add(self, element):
        self._mna.add_element(element)

    def dc(self):
        """Solve the linear DC operating point using the sources' DC values."""
        G, b_dc = self._mna.assemble()
        return Result(dc_solve(G, b_dc), self._mna.node_map)

    def ac(self, frequencies_hz, sources):
        """Sweep frequencies using {independent source name: complex phasor}.

        Omitted sources have zero AC excitation. Voltage phasors are in volts,
        current phasors in amps. The stored DC source values are independent.
        """
        G, M, _ = self._mna.assemble_system()
        nodes = self._mna.node_map
        b_ac = np.zeros(G.shape[0], dtype=complex)
        elements = {element.name: element for element in self._mna.elements}

        for name, phasor in sources.items():
            if name not in elements:
                raise ValueError(f"Unknown AC source {name!r}")
            source = elements[name]
            if isinstance(source, VoltageSource):
                b_ac[nodes.get_branch_index(name)] += phasor
            elif isinstance(source, CurrentSource):
                plus = nodes.get_index(source.node_plus)
                minus = nodes.get_index(source.node_minus)
                if plus is not None:
                    b_ac[plus] -= phasor
                if minus is not None:
                    b_ac[minus] += phasor
            else:
                raise ValueError(f"{name!r} must be an independent source")

        return Result(ac_sweep(G, M, b_ac, frequencies_hz), nodes)
