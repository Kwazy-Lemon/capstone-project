"""Map circuit node labels to Modified Nodal Analysis matrix indices."""


class NodeMap:
    """Assign one matrix index to each non-ground node."""

    def __init__(self):
        """Create an empty node-to-index mapping."""
        self.nodes = {}
        self.next_index = 0

    @staticmethod
    def _is_ground(node):
        """Return whether a label denotes the reference (ground) node."""
        return node.lower() in ("0", "gnd")

    def add_node(self, node):
        """Register a node and return its matrix index.

        Ground is the zero-volt reference, so it is not represented by an
        unknown in the MNA system. Repeated labels reuse their existing index.
        """
        if self._is_ground(node):
            return None

        if node not in self.nodes:
            self.nodes[node] = self.next_index
            self.next_index += 1

        return self.nodes[node]

    def get_index(self, node):
        """Return a node's index, or None for ground.

        Non-ground nodes must have been registered using add_node first.
        """
        if self._is_ground(node):
            return None

        return self.nodes[node]

    def num_nodes(self):
        """Return the number of registered non-ground nodes."""
        return self.next_index

    def __repr__(self):
        """Return a readable representation of the node-index mapping."""
        return str(self.nodes)
