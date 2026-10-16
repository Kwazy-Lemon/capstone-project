"""Ideal resistor model for DC Modified Nodal Analysis."""


class Resistor:
    """Represent a resistor connected between two circuit nodes."""

    def __init__(self, name, node1, node2, resistance):
        """Store the resistor label, terminal nodes, and resistance in ohms."""
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.resistance = resistance

    def stamp(self, node_map):
        """Return this resistor's contributions to the MNA matrix.

        Ohm's law gives the conductance g = 1 / R. For terminals a and b,
        the resistor contributes +g to each connected node's diagonal and -g
        to the off-diagonal entries when both terminals are non-ground.
        Ground entries are omitted because ground has no matrix index.

        Returns:
            List of (row, column, value) triples for the MNA matrix.
        """
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

        return entries
