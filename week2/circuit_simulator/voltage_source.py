class VoltageSource:
    def __init__(self, name, node_plus, node_minus, voltage):
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.voltage = voltage
        self.branch_index = None

    def assign_branch_index(self, index):
        """Assign the MNA index for the voltage-source current."""
        self.branch_index = index

    def stamp(self, node_map):
        """Return the MNA matrix and RHS contributions."""
        if self.branch_index is None:
            raise ValueError(
                f"Branch index has not been assigned to {self.name}."
            )

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)

        entries = []

        # KCL contribution of the voltage-source current
        if n_plus is not None:
            entries.append((n_plus, self.branch_index, 1.0))
            entries.append((self.branch_index, n_plus, 1.0))

        if n_minus is not None:
            entries.append((n_minus, self.branch_index, -1.0))
            entries.append((self.branch_index, n_minus, -1.0))

        # RHS contribution from the voltage constraint
        rhs_value = self.voltage

        return entries, (self.branch_index, rhs_value)
