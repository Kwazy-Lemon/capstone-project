"""Independent ideal voltage source model for DC Modified Nodal Analysis."""


class VoltageSource:
    """Represent an ideal source with positive and negative terminals."""

    def __init__(self, name, node_plus, node_minus, voltage):
        """Store the source label, terminal nodes, and voltage in volts."""
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.voltage = voltage
        # MNAAssembler assigns this index after all node-voltage indices exist.
        self.branch_index = None

    def assign_branch_index(self, index):
        """Set the solution-vector index for this source's current unknown."""
        self.branch_index = index

    def stamp(self, node_map):
        """Return matrix entries and the RHS entry for this voltage source.

        The ideal source imposes V(+) - V(-) = voltage. Its current cannot be
        represented by a conductance, so MNA adds a separate current unknown.
        Terminal coefficients appear in both the node KCL equations and the
        voltage-constraint equation.

        Returns:
            A pair (entries, rhs_entry). entries contains (row, column, value)
            triples; rhs_entry is the (row_index, value) for the voltage
            constraint.
        """
        if self.branch_index is None:
            raise ValueError(
                f"Branch index has not been assigned to {self.name}."
            )

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)

        entries = []

        # Couple the source-current unknown to the positive terminal's KCL
        # equation and to the source's voltage-constraint row.
        if n_plus is not None:
            entries.append((n_plus, self.branch_index, 1.0))
            entries.append((self.branch_index, n_plus, 1.0))

        # The negative terminal has the opposite sign in V(+) - V(-).
        if n_minus is not None:
            entries.append((n_minus, self.branch_index, -1.0))
            entries.append((self.branch_index, n_minus, -1.0))

        # Right-hand side of the constraint V(+) - V(-) = source voltage.
        rhs_entry = (self.branch_index, self.voltage)

        return entries, rhs_entry
