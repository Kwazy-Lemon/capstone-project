from math import isfinite
import numpy as np

class Diode:

    requires_branch_current = False

    def __init__(self, name, anode, cathode, Is=1e-14, n=1.0, Vt=0.02585):

        # raise error for diode parameters
        if not isfinite(Is) or Is <= 0:
            raise ValueError("Is must be finite and positive")
        if not isfinite(Vt) or Vt <= 0:
            raise ValueError("Vt must be finite and positive")
        if not isfinite(n) or n <= 0:
            raise ValueError("n must be finite and positive")

        self.name = name
        self.anode = anode
        self.cathode = cathode
        self.Is = Is
        self.n = n
        self.Vt = Vt

    @property
    def nodes(self):
        # return anode cathode terminals for node registration
        return (self.anode, self.cathode)


    def stamp(self, node_map):
        # return nothing as it will contribute nonlinear entries
        return [], [], []

    def nonlinear_contribution(self, node_map, x):
        anode_index = node_map.get_index(self.anode)
        cathode_index = node_map.get_index(self.cathode)

        # Ground node = 0 V condition
        V_A = 0.0 if anode_index is None else x[anode_index]
        V_C = 0.0 if cathode_index is None else x[cathode_index]

        # calculate diode current value
        Id = self.Is * np.expm1((V_A - V_C) / (self.n * self.Vt))

        # calculate differential conductance
        gd = (self.Is / (self.n * self.Vt)) * np.exp((V_A - V_C) / (self.n * self.Vt))

        # build f and jacobian entries
        f = []
        j = []

        if anode_index is not None:
            f.append((anode_index, Id))
            j.append((anode_index, anode_index, gd))

        if cathode_index is not None:
            f.append((cathode_index, -Id))
            j.append((cathode_index, cathode_index, gd))

        if anode_index is not None and cathode_index is not None:
            j.append((anode_index, cathode_index, -gd))
            j.append((cathode_index, anode_index, -gd))

        return f, j




