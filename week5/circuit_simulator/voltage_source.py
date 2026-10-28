"""Zhuoyao Luo"""

from math import isfinite


class VoltageSource:
    """represent a voltage source between two nodes"""

    # source current needs an extra unknown
    requires_branch_current = True

    def __init__(self, name, node_plus, node_minus, voltage):
        # voltage is measured in volts
        if not isfinite(voltage):
            raise ValueError("voltage must be finite")

        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.voltage = voltage

        # the assembler assigns this index before stamping
        self.branch_index = None

    @property
    def nodes(self):
        """return both terminals for node registration"""
        return (self.node_plus, self.node_minus)

    def assign_branch_index(self, index):
        """store the index of the source current"""
        self.branch_index = index

    def stamp(self, node_map):
        """return matrix and source contributions"""

        if self.branch_index is None:
            raise ValueError(
                f"branch index has not been assigned to {self.name!r}"
            )

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        k = self.branch_index

        g_entries = []

        # current leaves the positive terminal
        # the branch equation includes the positive terminal voltage
        if n_plus is not None:
            g_entries.append((n_plus, k, 1.0))
            g_entries.append((k, n_plus, 1.0))

        # current enters the negative terminal
        # the branch equation subtracts the negative terminal voltage
        if n_minus is not None:
            g_entries.append((n_minus, k, -1.0))
            g_entries.append((k, n_minus, -1.0))

        # set the voltage difference in the branch equation
        b_entries = [(k, self.voltage)]

        # an ideal voltage source has no derivative contribution
        return g_entries, [], b_entries
