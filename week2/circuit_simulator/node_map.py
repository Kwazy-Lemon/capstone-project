class NodeMap:
    def __init__(self):
        self.nodes = {}
        self.next_index = 0

    def add_node(self, node):
        """Add a node and return its matrix index."""
        if node.lower() in ("0", "gnd"):
            return None

        if node not in self.nodes:
            self.nodes[node] = self.next_index
            self.next_index += 1

        return self.nodes[node]

    def get_index(self, node):
        """Return the matrix index of a node."""
        if node.lower() in ("0", "gnd"):
            return None

        return self.nodes[node]

    def num_nodes(self):
        """Return the number of non-ground nodes."""
        return self.next_index

    def __repr__(self):
        return str(self.nodes)
