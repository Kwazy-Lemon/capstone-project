class Resistor:
    def __init__(self, name, node1, node2, resistance):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.resistance = resistance

    def stamp(self, node_map):
        """Return the MNA matrix contributions of this resistor."""
        conductance = 1.0 / self.resistance

        n1 = node_map.get_index(self.node1)
        n2 = node_map.get_index(self.node2)

        entries = []

        if n1 is not None:
            entries.append((n1, n1, conductance))

        if n2 is not None:
            entries.append((n2, n2, conductance))

        if n1 is not None and n2 is not None:
            entries.append((n1, n2, -conductance))
            entries.append((n2, n1, -conductance))

        return entries
