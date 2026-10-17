"""independent ideal voltage source model for DC Modified Nodal Analysis"""


class VoltageSource:
    """Represent an ideal source with positive and negative terminals"""

    requires_branch_current = True

    def __init__(self, name, node_plus, node_minus, voltage):
        """store the source label, terminal nodes, and voltage in volts"""
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.voltage = voltage
        # MNAAssembler assigns this index after all node-voltage indices exist.
        self.branch_index = None

    @property
    def nodes(self):
        return (self.node_plus, self.node_minus)

    def assign_branch_index(self, index):
        """set the solution-vector index for this source's current unknown"""
        self.branch_index = index

    def stamp(self, node_map):
        """return (G entries, M entries, RHS entries) for V(+) - V(-) = V"""
        if self.branch_index is None:
            raise ValueError(
                f"Branch index has not been assigned to {self.name}."
            )

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)

        entries = []

        # Couple branch current and terminal voltage.
        if n_plus is not None:
            entries.append((n_plus, self.branch_index, 1.0))
            entries.append((self.branch_index, n_plus, 1.0))

        # The negative terminal has the opposite sign in V(+) - V(-).
        if n_minus is not None:
            entries.append((n_minus, self.branch_index, -1.0))
            entries.append((self.branch_index, n_minus, -1.0))

        # Right-hand side of the constraint V(+) - V(-) = source voltage.
        rhs_entry = (self.branch_index, self.voltage)

        return entries, [], [rhs_entry]
