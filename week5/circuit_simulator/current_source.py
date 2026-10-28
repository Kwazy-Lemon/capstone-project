"""Zhuoyao Luo"""

from math import isfinite


class CurrentSource:
    """represent a current source between two nodes"""

    # source current is known and needs no extra unknown
    requires_branch_current = False

    def __init__(self, name, node_plus, node_minus, current):
        # current is measured in amperes
        if not isfinite(current):
            raise ValueError("current must be finite")

        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.current = current

    @property
    def nodes(self):
        """return both terminals for node registration"""
        return (self.node_plus, self.node_minus)

    def stamp(self, node_map):
        """return matrix and source contributions"""

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)

        b_entries = []

        # current leaves the positive terminal
        if n_plus is not None:
            b_entries.append((n_plus, -self.current))

        # current enters the negative terminal
        if n_minus is not None:
            b_entries.append((n_minus, self.current))

        # a current source only contributes to the source vector
        return [], [], b_entries
