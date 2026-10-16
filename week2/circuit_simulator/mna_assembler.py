"""Assemble the Modified Nodal Analysis system for the supported circuit elements.

The assembled linear system is A @ x = rhs. Its unknown vector contains the
voltages of non-ground nodes followed by the currents through ideal voltage
sources.
"""

import numpy as np
from scipy.sparse import coo_matrix

from node_map import NodeMap
from voltage_source import VoltageSource


class MNAAssembler:
    """Collect circuit elements and construct a sparse MNA matrix."""

    def __init__(self):
        """Initialize node tracking and element collections."""
        self.node_map = NodeMap()
        self.elements = []
        self.voltage_sources = []

    def add_element(self, element):
        """Register an element and add its terminal nodes to the node map.

        Resistors use node1/node2 attributes; voltage sources use
        node_plus/node_minus. Ground is recognized by NodeMap and is not
        assigned a voltage unknown.
        """
        self.elements.append(element)

        # Register terminals before assigning branch-current indices.
        if hasattr(element, "node1"):
            self.node_map.add_node(element.node1)
            self.node_map.add_node(element.node2)
        elif hasattr(element, "node_plus"):
            self.node_map.add_node(element.node_plus)
            self.node_map.add_node(element.node_minus)

        # Each ideal voltage source introduces one extra current unknown.
        if isinstance(element, VoltageSource):
            self.voltage_sources.append(element)

    def assign_branch_indices(self):
        """Assign matrix indices to voltage-source current unknowns.

        Branch currents are placed after the non-ground node voltages in the
        solution vector.
        """
        next_index = self.node_map.num_nodes()

        for voltage_source in self.voltage_sources:
            voltage_source.assign_branch_index(next_index)
            next_index += 1

    def assemble(self):
        """Return the sparse MNA matrix and right-hand-side vector.

        Each element's stamp method returns its contributions. COO format is
        convenient while collecting entries; conversion to CSC produces the
        sparse matrix returned to the caller.

        Returns:
            tuple: (matrix, rhs), representing matrix @ x = rhs.
        """
        # Assign source-current indices after all nodes have been registered.
        self.assign_branch_indices()

        # Unknowns = non-ground node voltages + voltage-source currents.
        n = self.node_map.num_nodes() + len(self.voltage_sources)

        # COO stores each contribution as a row, column, and value triple.
        # Repeated entries at the same position are combined by sparse format
        # conversion, which lets separate elements contribute to one entry.
        rows = []
        cols = []
        values = []
        rhs = np.zeros(n)

        for element in self.elements:
            contribution = element.stamp(self.node_map)

            if isinstance(element, VoltageSource):
                # A voltage source contributes matrix entries and one voltage
                # constraint value in the right-hand-side vector.
                entries, rhs_entry = contribution
                for row, col, value in entries:
                    rows.append(row)
                    cols.append(col)
                    values.append(value)

                rhs_index, rhs_value = rhs_entry
                rhs[rhs_index] += rhs_value
            else:
                # Resistors contribute conductance-matrix entries only.
                for row, col, value in contribution:
                    rows.append(row)
                    cols.append(col)
                    values.append(value)

        matrix = coo_matrix(
            (values, (rows, cols)),
            shape=(n, n),
        ).tocsc()

        return matrix, rhs
