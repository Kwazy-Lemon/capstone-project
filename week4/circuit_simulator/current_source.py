"""Independent ideal current source for the shared component interface."""

from math import isfinite


class CurrentSource:
    """A prescribed real current, in amperes, from node_plus to node_minus.

    A negative value reverses the actual current direction. Zero is allowed.
    As with the existing VoltageSource, the stored value is constant. AC
    phasors and transient waveforms can be supplied through the solver's RHS.
    """

    # Current is already known, so no extra branch-current unknown is needed.
    requires_branch_current = False

    def __init__(self, name, node_plus, node_minus, current):
        if not isfinite(current):
            raise ValueError("Source current must be finite (amperes).")
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.current = current

    @property
    def nodes(self):
        """Register both terminals even when this source is added first."""
        return (self.node_plus, self.node_minus)

    def stamp(self, node_map):
        """Return RHS injections; this ideal source changes neither G nor M.

        Node equations use 'passive current leaving the node = injected
        source current'. A positive source current leaves node_plus, so it
        contributes -I there. It enters node_minus, contributing +I there.
        """
        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        rhs_entries = []

        if n_plus is not None:
            rhs_entries.append((n_plus, -self.current))
        if n_minus is not None:
            rhs_entries.append((n_minus, self.current))

        # Ground needs no equation. If the two labels refer to the same node,
        # the assembler sums -I and +I and they correctly cancel.
        return [], [], rhs_entries
