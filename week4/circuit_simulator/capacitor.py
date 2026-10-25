"""Ideal capacitor using the existing G*x + M*dx/dt = b interface.

Terminal voltage is V(node_plus) - V(node_minus). Current is positive from
node_plus to node_minus, so the passive sign convention gives i = C*dv/dt.
"""

from math import isfinite


class Capacitor:
    """A constant, positive capacitance measured in farads (F)."""

    # Capacitor current can be written using voltage derivatives. Unlike an
    # inductor, this component does not need a separate branch-current unknown.
    requires_branch_current = False

    def __init__(self, name, node_plus, node_minus, capacitance):
        if not isfinite(capacitance) or capacitance <= 0:
            raise ValueError("Capacitance must be finite and positive (farads).")
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.capacitance = capacitance

    @property
    def nodes(self):
        """Expose terminal labels for registration by MNAAssembler."""
        return (self.node_plus, self.node_minus)

    def stamp(self, node_map):
        """Return derivative-matrix entries from i = C*d(Vplus - Vminus)/dt.

        The capacitor leaves G and b unchanged. At DC the derivative is zero,
        so this stamp naturally describes an open circuit. The transient
        solver receives initial capacitor voltage through its node-voltage x0.
        """
        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        c = self.capacitance
        m_entries = []

        # KCL counts currents leaving each node. At the positive terminal,
        # i_C = +C*dVplus/dt - C*dVminus/dt, giving the first diagonal +C.
        if n_plus is not None:
            m_entries.append((n_plus, n_plus, c))

        # At the negative terminal, -i_C = +C*dVminus/dt - C*dVplus/dt.
        # Therefore the second diagonal is positive too.
        if n_minus is not None:
            m_entries.append((n_minus, n_minus, c))

        # If both nodes have unknown voltages, their derivatives are coupled.
        # Ground is fixed at zero, so its derivative and entries are omitted.
        if n_plus is not None and n_minus is not None:
            m_entries.append((n_plus, n_minus, -c))
            m_entries.append((n_minus, n_plus, -c))

        return [], m_entries, []
