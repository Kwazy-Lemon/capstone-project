import numpy as np
from scipy.sparse import coo_matrix
from node_map import NodeMap
from voltage_source import VoltageSource


class MNAAssembler:
    def __init__(self):
        self.node_map = NodeMap()
        self.elements = []
        self.voltage_sources = []

    def add_element(self, element):
        """Add a circuit element to the assembler."""
        self.elements.append(element)

        # Register the element's nodes
        if hasattr(element, "node1"):
            self.node_map.add_node(element.node1)
            self.node_map.add_node(element.node2)

        elif hasattr(element, "node_plus"):
            self.node_map.add_node(element.node_plus)
            self.node_map.add_node(element.node_minus)

        # Keep track of voltage sources
        if isinstance(element, VoltageSource):
            self.voltage_sources.append(element)

    def assign_branch_indices(self):
        """Assign MNA indices to voltage-source currents."""
        next_index = self.node_map.num_nodes()

        for voltage_source in self.voltage_sources:
            voltage_source.assign_branch_index(next_index)
            next_index += 1

    def assemble(self):
        """Assemble the MNA matrix and RHS vector."""
        self.assign_branch_indices()

        n = self.node_map.num_nodes() + len(self.voltage_sources)

        rows = []
        cols = []
        values = []

        rhs = np.zeros(n)

        for element in self.elements:
            contribution = element.stamp(self.node_map)

            if isinstance(element, VoltageSource):
                entries, rhs_entry = contribution

                for row, col, value in entries:
                    rows.append(row)
                    cols.append(col)
                    values.append(value)

                rhs_index, rhs_value = rhs_entry
                rhs[rhs_index] += rhs_value

            else:
                entries = contribution

                for row, col, value in entries:
                    rows.append(row)
                    cols.append(col)
                    values.append(value)

        matrix = coo_matrix(
            (values, (rows, cols)),
            shape=(n, n)
        ).tocsc()

        return matrix, rhs
