"""Zhuoyao Luo"""

from math import isfinite


class Resistor:
    """represent a resistor between two nodes"""

    # resistor current does not need an extra unknown
    requires_branch_current = False

    def __init__(self, name, node1, node2, resistance):
        # resistance is measured in ohms
        if not isfinite(resistance) or resistance <= 0:
            raise ValueError("resistance must be finite and positive")

        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.resistance = resistance

    @property
    def nodes(self):
        """return both terminals for node registration"""
        return (self.node1, self.node2)

    def stamp(self, node_map):
        """return matrix and source contributions"""

        conductance = 1.0 / self.resistance

        n1 = node_map.get_index(self.node1)
        n2 = node_map.get_index(self.node2)

        g_entries = []

        # add diagonal entries for non-ground nodes
        if n1 is not None:
            g_entries.append((n1, n1, conductance))

        if n2 is not None:
            g_entries.append((n2, n2, conductance))

        # add coupling entries when neither terminal is ground
        if n1 is not None and n2 is not None:
            g_entries.append((n1, n2, -conductance))
            g_entries.append((n2, n1, -conductance))

        # a resistor has no derivative or source contribution
        return g_entries, [], []
