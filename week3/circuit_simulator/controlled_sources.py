"""Four ideal linear controlled sources with constant, real gains.

Voltage differences are always V(+) - V(-). Output and control branch
currents are positive from the corresponding positive to negative terminal.
Each class writes its equation directly; no additional class hierarchy is
needed. All four devices are algebraic, so their M and RHS entries are empty.
"""

from math import isfinite


class VCVS:
    """Voltage-controlled voltage source: Vout = gain * Vcontrol.

    gain is a dimensionless voltage ratio (V/V), including zero or negative
    values. Control terminals sense voltage without drawing any input current.
    """

    requires_branch_current = True

    def __init__(self, name, node_plus, node_minus,
                 control_plus, control_minus, gain):
        if not isfinite(gain):
            raise ValueError("VCVS gain must be finite.")
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.control_plus = control_plus
        self.control_minus = control_minus
        self.gain = gain
        self.branch_index = None

    @property
    def nodes(self):
        # Control nodes must also be registered, even if this device is added
        # before the components that drive its control voltage.
        return (self.node_plus, self.node_minus,
                self.control_plus, self.control_minus)

    def assign_branch_index(self, index):
        self.branch_index = index

    def stamp(self, node_map):
        """Write V(+) - V(-) - gain*(Vcontrol(+) - Vcontrol(-)) = 0."""
        if self.branch_index is None:
            raise ValueError(f"Branch index has not been assigned to {self.name}.")
        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        c_plus = node_map.get_index(self.control_plus)
        c_minus = node_map.get_index(self.control_minus)
        k = self.branch_index
        entries = []

        # The unknown source current couples to output KCL, just as for an
        # independent voltage source. The new row imposes its voltage law.
        if n_plus is not None:
            entries.append((n_plus, k, 1.0))
            entries.append((k, n_plus, 1.0))
        if n_minus is not None:
            entries.append((n_minus, k, -1.0))
            entries.append((k, n_minus, -1.0))

        # Expand -gain*(Vcontrol(+) - Vcontrol(-)). Only the constraint
        # row receives these terms: an ideal voltage input draws no current.
        if c_plus is not None:
            entries.append((k, c_plus, -self.gain))
        if c_minus is not None:
            entries.append((k, c_minus, self.gain))
        return entries, [], []


class VCCS:
    """Voltage-controlled current source: Iout = gm * Vcontrol.

    transconductance (gm) is in siemens (A/V). The specified current flows
    from node_plus to node_minus; a negative value reverses that direction.
    """

    # Iout is already expressed using node voltages, so no extra unknown is
    # necessary. This differs from both kinds of voltage-output source.
    requires_branch_current = False

    def __init__(self, name, node_plus, node_minus,
                 control_plus, control_minus, transconductance):
        if not isfinite(transconductance):
            raise ValueError("VCCS transconductance must be finite.")
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.control_plus = control_plus
        self.control_minus = control_minus
        self.transconductance = transconductance

    @property
    def nodes(self):
        return (self.node_plus, self.node_minus,
                self.control_plus, self.control_minus)

    def stamp(self, node_map):
        """Add +Iout to positive-terminal KCL and -Iout to negative KCL."""
        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        c_plus = node_map.get_index(self.control_plus)
        c_minus = node_map.get_index(self.control_minus)
        gm = self.transconductance
        entries = []

        # Positive output node: +gm*Vcontrol(+) - gm*Vcontrol(-).
        if n_plus is not None:
            if c_plus is not None:
                entries.append((n_plus, c_plus, gm))
            if c_minus is not None:
                entries.append((n_plus, c_minus, -gm))

        # The same current enters the negative output node, hence opposite
        # signs. Omitting ground rows/columns works for all four terminals.
        if n_minus is not None:
            if c_plus is not None:
                entries.append((n_minus, c_plus, -gm))
            if c_minus is not None:
                entries.append((n_minus, c_minus, gm))
        return entries, [], []


class CCCS:
    """Current-controlled current source: Iout = gain * Icontrol.

    gain is a dimensionless current ratio (A/A). control_branch is the exact
    name of a device with an allocated current unknown, e.g. a voltage source
    or inductor. To sense a resistor current, insert a 0 V source in series.
    """

    requires_branch_current = False

    def __init__(self, name, node_plus, node_minus, control_branch, gain):
        if not isfinite(gain):
            raise ValueError("CCCS gain must be finite.")
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.control_branch = control_branch
        self.gain = gain

    @property
    def nodes(self):
        # A branch name is not a node label. Its device registers its own nodes.
        return (self.node_plus, self.node_minus)

    def stamp(self, node_map):
        """Couple output KCL to the named control-current column."""
        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)

        # All branch indices exist before stamping, so forward references are
        # valid. A missing branch gives the existing NodeMap's clear ValueError.
        control = node_map.get_branch_index(self.control_branch)
        entries = []
        if n_plus is not None:
            entries.append((n_plus, control, self.gain))
        if n_minus is not None:
            entries.append((n_minus, control, -self.gain))
        return entries, [], []


class CCVS:
    """Current-controlled voltage source: Vout = rm * Icontrol.

    transresistance (rm) is in ohms (V/A). The output needs a new branch
    current; control_branch refers to a separately allocated current unknown.
    """

    requires_branch_current = True

    def __init__(self, name, node_plus, node_minus,
                 control_branch, transresistance):
        if not isfinite(transresistance):
            raise ValueError("CCVS transresistance must be finite.")
        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.control_branch = control_branch
        self.transresistance = transresistance
        self.branch_index = None

    @property
    def nodes(self):
        return (self.node_plus, self.node_minus)

    def assign_branch_index(self, index):
        self.branch_index = index

    def stamp(self, node_map):
        """Write V(+) - V(-) - rm*Icontrol = 0, plus output KCL terms."""
        if self.branch_index is None:
            raise ValueError(f"Branch index has not been assigned to {self.name}.")
        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        control = node_map.get_branch_index(self.control_branch)
        k = self.branch_index
        entries = []

        # These entries are the same output-current/voltage coupling used by
        # VCVS. The control law below replaces an independent source voltage.
        if n_plus is not None:
            entries.append((n_plus, k, 1.0))
            entries.append((k, n_plus, 1.0))
        if n_minus is not None:
            entries.append((n_minus, k, -1.0))
            entries.append((k, n_minus, -1.0))

        # The control current is an unknown in x, so its coefficient belongs
        # in G, not in b. Keeping it on the left gives the minus sign.
        entries.append((k, control, -self.transresistance))
        return entries, [], []
