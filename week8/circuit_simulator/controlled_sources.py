"""Zhuoyao Luo"""

from math import isfinite


class VCVS:
    """represent a voltage controlled voltage source"""

    # output voltage requires a branch current unknown
    requires_branch_current = True

    def __init__(
        self, name, node_plus, node_minus,
        control_plus, control_minus, gain
    ):
        if not isfinite(gain):
            raise ValueError("gain must be finite")

        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.control_plus = control_plus
        self.control_minus = control_minus
        self.gain = gain
        self.branch_index = None

    @property
    def nodes(self):
        """return output and control nodes"""
        return (
            self.node_plus,
            self.node_minus,
            self.control_plus,
            self.control_minus,
        )

    def assign_branch_index(self, index):
        """store the output current index"""
        self.branch_index = index

    def stamp(self, node_map):
        """return the voltage source contributions"""

        if self.branch_index is None:
            raise ValueError(
                f"branch index has not been assigned to {self.name!r}"
            )

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        c_plus = node_map.get_index(self.control_plus)
        c_minus = node_map.get_index(self.control_minus)
        k = self.branch_index

        g_entries = []

        # add output current and voltage terms
        if n_plus is not None:
            g_entries.append((n_plus, k, 1.0))
            g_entries.append((k, n_plus, 1.0))

        if n_minus is not None:
            g_entries.append((n_minus, k, -1.0))
            g_entries.append((k, n_minus, -1.0))

        # subtract the controlled voltage from the branch equation
        if c_plus is not None:
            g_entries.append((k, c_plus, -self.gain))

        if c_minus is not None:
            g_entries.append((k, c_minus, self.gain))

        return g_entries, [], []


class VCCS:
    """represent a voltage controlled current source"""

    # output current is expressed using node voltages
    requires_branch_current = False

    def __init__(
        self, name, node_plus, node_minus,
        control_plus, control_minus, transconductance
    ):
        if not isfinite(transconductance):
            raise ValueError("transconductance must be finite")

        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.control_plus = control_plus
        self.control_minus = control_minus
        self.transconductance = transconductance

    @property
    def nodes(self):
        """return output and control nodes"""
        return (
            self.node_plus,
            self.node_minus,
            self.control_plus,
            self.control_minus,
        )

    def stamp(self, node_map):
        """return the current source contributions"""

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        c_plus = node_map.get_index(self.control_plus)
        c_minus = node_map.get_index(self.control_minus)
        gm = self.transconductance

        g_entries = []

        # add output current to the positive node equation
        if n_plus is not None:
            if c_plus is not None:
                g_entries.append((n_plus, c_plus, gm))
            if c_minus is not None:
                g_entries.append((n_plus, c_minus, -gm))

        # subtract output current from the negative node equation
        if n_minus is not None:
            if c_plus is not None:
                g_entries.append((n_minus, c_plus, -gm))
            if c_minus is not None:
                g_entries.append((n_minus, c_minus, gm))

        return g_entries, [], []


class CCCS:
    """represent a current controlled current source"""

    # output current uses an existing branch current unknown
    requires_branch_current = False

    def __init__(self, name, node_plus, node_minus, control_branch, gain):
        if not isfinite(gain):
            raise ValueError("gain must be finite")

        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.control_branch = control_branch
        self.gain = gain

    @property
    def nodes(self):
        """return the output nodes"""
        return (self.node_plus, self.node_minus)

    def stamp(self, node_map):
        """return the current source contributions"""

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        control = node_map.get_branch_index(self.control_branch)

        g_entries = []

        # add controlled current to the output node equations
        if n_plus is not None:
            g_entries.append((n_plus, control, self.gain))

        if n_minus is not None:
            g_entries.append((n_minus, control, -self.gain))

        return g_entries, [], []


class CCVS:
    """represent a current controlled voltage source"""

    # output voltage requires a branch current unknown
    requires_branch_current = True

    def __init__(
        self, name, node_plus, node_minus,
        control_branch, transresistance
    ):
        if not isfinite(transresistance):
            raise ValueError("transresistance must be finite")

        self.name = name
        self.node_plus = node_plus
        self.node_minus = node_minus
        self.control_branch = control_branch
        self.transresistance = transresistance
        self.branch_index = None

    @property
    def nodes(self):
        """return the output nodes"""
        return (self.node_plus, self.node_minus)

    def assign_branch_index(self, index):
        """store the output current index"""
        self.branch_index = index

    def stamp(self, node_map):
        """return the voltage source contributions"""

        if self.branch_index is None:
            raise ValueError(
                f"branch index has not been assigned to {self.name!r}"
            )

        n_plus = node_map.get_index(self.node_plus)
        n_minus = node_map.get_index(self.node_minus)
        control = node_map.get_branch_index(self.control_branch)
        k = self.branch_index

        g_entries = []

        # add output current and voltage terms
        if n_plus is not None:
            g_entries.append((n_plus, k, 1.0))
            g_entries.append((k, n_plus, 1.0))

        if n_minus is not None:
            g_entries.append((n_minus, k, -1.0))
            g_entries.append((k, n_minus, -1.0))

        # subtract the controlled current term from the branch equation
        g_entries.append((k, control, -self.transresistance))

        return g_entries, [], []
