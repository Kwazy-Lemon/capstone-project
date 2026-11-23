"""Zhuoyao Luo"""

from math import isfinite


class Inductor:
    """represent an inductor between two nodes"""

    # inductor current needs an extra unknown
    requires_branch_current = True

    def __init__(self, name, node_plus, node_minus, inductance):
        # inductance is measured in henries
        if not isfinite(inductance) or inductance <= 0:
            raise ValueError("inductance must be finite and positive")

        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.inductance = inductance

        # the assembler assigns the branch index before stamping
        self.branch_index = None

    @property
    def nodes(self):
        """return both terminals for node registration"""
        return (self.node_plus, self.node_minus)

    def assign_branch_index(self, index):
        """store the index of the inductor current"""
        self.branch_index = index

    def stamp(self, node_map):
        """return the inductor contributions"""

        if self.branch_index is None:
            raise ValueError(
                f"branch index has not been assigned to {self.name!r}"
            )

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        k = self.branch_index

        g_entries = []

        # add current and voltage terms for the positive terminal
        if n_plus is not None:
            g_entries.append((n_plus, k, 1.0))
            g_entries.append((k, n_plus, 1.0))

        # add current and voltage terms for the negative terminal
        if n_minus is not None:
            g_entries.append((n_minus, k, -1.0))
            g_entries.append((k, n_minus, -1.0))

        # add the current derivative term
        m_entries = [(k, k, -self.inductance)]

        # an inductor has no independent source contribution
        return g_entries, m_entries, []
