"""map circuit node labels to Modified Nodal Analysis matrix indices"""


class NodeMap:
    """assign one matrix index to each non-ground node"""

    def __init__(self):
        """create an empty node-to-index mapping"""
        self.nodes = {}
        self.branches = {}
        self.next_index = 0

    @staticmethod
    def _is_ground(node):
        """return whether a label denotes the reference (ground) node"""
        return node.lower() in ("0", "gnd")

    def add_node(self, node):
        """register a node and return its matrix index.

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
        """return a node's index, or None for ground

        Non-ground nodes must have been registered using add_node first
        """
        if self._is_ground(node):
            return None

        return self.nodes[node]

    def num_nodes(self):
        """Return the number of registered non-ground nodes."""
        return self.next_index

    def get_branch_index(self, name):
        """Look up a named current unknown, including a control branch."""
        try:
            return self.branches[name]
        except KeyError as error:
            raise ValueError(f"No branch-current unknown named {name!r}.") from error

    def __repr__(self):
        """Return a readable representation of the node-index mapping."""
        return str(self.nodes)
