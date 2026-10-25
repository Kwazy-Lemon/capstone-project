"""ideal resistor model for DC Modified Nodal Analysis"""


class Resistor:
    """represent a resistor connected between two circuit nodes"""

    requires_branch_current = False

    def __init__(self, name, node1, node2, resistance):
        """store the resistor label, terminal nodes, and resistance in ohms"""
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.resistance = resistance

    @property
    def nodes(self):
        return (self.node1, self.node2)

    def stamp(self, node_map):
        """return (G entries, M entries, RHS entries)"""
        conductance = 1.0 / self.resistance

        n1 = node_map.get_index(self.node1)
        n2 = node_map.get_index(self.node2)

        entries = []

        # Diagonal terms for each non-ground terminal.
        if n1 is not None:
            entries.append((n1, n1, conductance))

        if n2 is not None:
            entries.append((n2, n2, conductance))

        # Off-diagonal coupling exists only when both terminals are unknowns.
        if n1 is not None and n2 is not None:
            entries.append((n1, n2, -conductance))
            entries.append((n2, n1, -conductance))

        return entries, [], []
