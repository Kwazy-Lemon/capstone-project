"""Zhuoyao Luo"""

import numpy as np
from scipy.sparse import coo_matrix

from node_map import NodeMap


class MNAAssembler:
    """collect components and assemble circuit equations"""

    def __init__(self):
        self.node_map = NodeMap()
        self.elements = []
        self.branch_elements = []

    def add_element(self, element):
        """add a component and register its nodes"""

        # component names must be unique
        for existing in self.elements:
            if existing.name == element.name:
                raise ValueError(
                    f"duplicate component name {element.name!r}"
                )

        self.elements.append(element)

        # register all nodes used by this component
        for node in element.nodes:
            self.node_map.add_node(node)

        # remember components that need a current unknown
        if element.requires_branch_current:
            self.branch_elements.append(element)

    def assign_branch_indices(self):
        """assign current indices after all node voltage indices"""

        self.node_map.branches.clear()
        next_index = self.node_map.num_nodes()

        # assign every branch before any component is stamped
        for element in self.branch_elements:
            element.assign_branch_index(next_index)
            self.node_map.branches[element.name] = next_index
            next_index += 1

    def assemble(self):
        """return the matrix and source vector for dc analysis"""

        g_matrix, _, b = self.assemble_system()
        return g_matrix, b

    def assemble_system(self):
        """return conductance matrix derivative matrix and source vector"""

        self.assign_branch_indices()

        # count node voltages and extra branch currents
        size = (
            self.node_map.num_nodes()
            + len(self.branch_elements)
        )

        g_entries = []
        m_entries = []
        b_entries = []

        # collect contributions from every component
        for element in self.elements:
            g, m, b = element.stamp(self.node_map)

            g_entries.extend(g)
            m_entries.extend(m)
            b_entries.extend(b)

        g_matrix = self._build_matrix(g_entries, size)
        m_matrix = self._build_matrix(m_entries, size)

        # add source contributions that share the same row
        b = np.zeros(size)

        for row, value in b_entries:
            b[row] += value

        return g_matrix, m_matrix, b

    @staticmethod
    def _build_matrix(entries, size):
        """convert row column value entries into a sparse matrix"""

        if not entries:
            return coo_matrix(
                (size, size), dtype=float
            ).tocsc()

        rows, columns, values = zip(*entries)

        # conversion adds contributions at matching positions
        return coo_matrix(
            (values, (rows, columns)),
            shape=(size, size),
            dtype=float,
        ).tocsc()
