"""Zhuoyao Luo"""

from math import isfinite


class Capacitor:
    """represent a capacitor between two nodes"""

    # capacitor current is expressed using voltage derivatives
    requires_branch_current = False

    def __init__(self, name, node_plus, node_minus, capacitance):
        # capacitance is measured in farads
        if not isfinite(capacitance) or capacitance <= 0:
            raise ValueError("capacitance must be finite and positive")

        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.capacitance = capacitance

    @property
    def nodes(self):
        """return both terminals for node registration"""
        return (self.node_plus, self.node_minus)

    def stamp(self, node_map):
        """return matrix and source contributions"""

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        c = self.capacitance

        m_entries = []

        # add diagonal entries for non-ground nodes
        if n_plus is not None:
            m_entries.append((n_plus, n_plus, c))

        if n_minus is not None:
            m_entries.append((n_minus, n_minus, c))

        # add coupling entries when neither terminal is ground
        if n_plus is not None and n_minus is not None:
            m_entries.append((n_plus, n_minus, -c))
            m_entries.append((n_minus, n_plus, -c))

        # a capacitor only contributes to the derivative matrix
        return [], m_entries, []
