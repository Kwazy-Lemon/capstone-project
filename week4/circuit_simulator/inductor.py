"""Ideal inductor using the shared Week 3 MNA component interface.

The branch current is positive from node_plus to node_minus. With this
reference direction, the passive sign convention gives v = L * di/dt.
"""

from math import isfinite


class Inductor:
    """An ideal, uncoupled inductor; inductance is measured in henries (H)."""

    # Unlike a resistor, an inductor needs its own current in the unknown vector.
    requires_branch_current = True

    def __init__(self, name, node_plus, node_minus, inductance):
        if not isfinite(inductance) or inductance <= 0:
            raise ValueError("Inductance must be finite and positive (henries).")
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.inductance = inductance
        self.branch_index = None

    @property
    def nodes(self):
        """Expose both terminals so the assembler can register their voltages."""
        return (self.node_plus, self.node_minus)

    def assign_branch_index(self, index):
        """The assembler allocates current indices after all node voltages."""
        self.branch_index = index

    def stamp(self, node_map):
        """Return contributions to G*x + M*dx/dt = b.

        The extra branch equation is V(+) - V(-) - L*dI_L/dt = 0.
        Thus G contains the terminal voltage terms, while M[k, k] = -L.
        The negative M entry follows from moving L*dI_L/dt to the left;
        it does not mean the physical inductance is negative.
        """
        if self.branch_index is None:
            raise ValueError(f"Branch index has not been assigned to {self.name}.")

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        k = self.branch_index
        g_entries = []

        # In node KCL equations we count currents leaving the node as positive.
        # I_L leaves the positive terminal, so its coefficient is +1 there.
        # The matching branch-row entry contributes +V(+) to the voltage law.
        if n_plus is not None:
            g_entries.append((n_plus, k, 1.0))
            g_entries.append((k, n_plus, 1.0))

        # At the negative terminal I_L enters the node; both signs are reversed.
        # Ground has no unknown, so entries involving its index are omitted.
        if n_minus is not None:
            g_entries.append((n_minus, k, -1.0))
            g_entries.append((k, n_minus, -1.0))

        m_entries = [(k, k, -self.inductance)]

        # This law has no independent forcing term. Initial current belongs to
        # the transient solver's x0, not to a permanent source contribution.
        # DC ignores M and obtains V(+) = V(-), an ideal short circuit.
        return g_entries, m_entries, []
