"""Assemble G*x + M*dx/dt = rhs from component stamps."""

import numpy as np
from scipy.sparse import coo_matrix

from node_map import NodeMap


class MNAAssembler:
    """collect elements and build sparse MNA matrices"""

    def __init__(self):
        self.node_map = NodeMap()
        self.elements = []
        self.branch_elements = []

    def add_element(self, element):
        """register output/control nodes and any extra branch current"""
        if any(existing.name == element.name for existing in self.elements):
            raise ValueError(f"Duplicate component name: {element.name!r}")
        self.elements.append(element)
        for node in element.nodes:
            self.node_map.add_node(node)
        if element.requires_branch_current:
            self.branch_elements.append(element)

    def assign_branch_indices(self):
        """allocate every branch before stamping to allow forward references"""
        self.node_map.branches.clear()
        next_index = self.node_map.num_nodes()
        for element in self.branch_elements:
            element.assign_branch_index(next_index)
            self.node_map.branches[element.name] = next_index
            next_index += 1

    def assemble(self):
        """keep the Week 2 DC interface: return (G, rhs)"""
        matrix, _, rhs = self.assemble_system()
        return matrix, rhs

    def assemble_system(self):
        """return (G, M, rhs) for DC, AC or transient analysis"""
        self.assign_branch_indices()
        n = self.node_map.num_nodes() + len(self.branch_elements)
        g_entries, m_entries, rhs_entries = [], [], []

        for element in self.elements:
            g, m, rhs = element.stamp(self.node_map)
            g_entries.extend(g)
            m_entries.extend(m)
            rhs_entries.extend(rhs)

        def sparse(entries):
            if not entries:
                return coo_matrix((n, n), dtype=float).tocsc()
            rows, cols, values = zip(*entries)
            # CSC conversion sums entries shared by multiple elements.
            return coo_matrix(
                (np.asarray(values), (rows, cols)), shape=(n, n)
            ).tocsc()

        rhs = np.zeros(n)
        for row, value in rhs_entries:
            rhs[row] += value
        return sparse(g_entries), sparse(m_entries), rhs

