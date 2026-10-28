"""Zhuoyao Luo"""


class NodeMap:
    """store indices used by the circuit matrix"""

    def __init__(self):
        # store non-ground node indices
        self.nodes = {}

        # the assembler fills in branch current indices
        self.branches = {}

        # index for the next new node
        self.next_index = 0

    @staticmethod
    def _is_ground(node):
        # accept ground labels regardless of letter case
        return node.lower() in ("0", "gnd")

    def add_node(self, node):
        """register a node and return its index"""

        # ground voltage is known and needs no matrix index
        if self._is_ground(node):
            return None

        # assign an index only when the node is new
        if node not in self.nodes:
            self.nodes[node] = self.next_index
            self.next_index += 1

        return self.nodes[node]

    def get_index(self, node):
        """return the index of a registered node"""

        if self._is_ground(node):
            return None

        if node not in self.nodes:
            raise ValueError(f"node {node!r} has not been registered")

        return self.nodes[node]

    def num_nodes(self):
        """return the number of non-ground nodes"""
        return self.next_index

    def get_branch_index(self, name):
        """return the index of a branch current"""

        if name not in self.branches:
            raise ValueError(f"no branch current unknown named {name!r}")

        return self.branches[name]

    def __repr__(self):
        """show the node labels and their indices"""
        return str(self.nodes)
