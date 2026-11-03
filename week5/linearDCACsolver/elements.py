from dataclasses import dataclass 
"""this file contains the basic dc elements like resistor, capacitor, inductor, 
VCVS VCCS CCVS CCCS."""

@dataclass
class Resistor:
    name: str
    node1: str
    node2: str
    value: float
@dataclass
class Capacitor:
    name: str
    node1: str
    node2: str
    value: float
@dataclass
class Inductor:
    name: str
    node1: str
    node2: str
    value: float

"""yes they all have 2 nodes, 1 value, and a name."""

@dataclass
class VoltageSource:
    name: str
    node1: str
    node2: str
    value: float

@dataclass
class CurrentSource:
    name: str
    node1: str
    node2: str
    value: float
"""normal voltage and current sources have 2 nodes, 1 value, and a name."""



"""
VCVS vout1 - vout2 = gain * (v_ctrl1 - v_ctrl2)

VCCS i_out = gain * (v_ctrl1 - v_ctrl2)

CCVS v_out = gain * i_ctrl

CCCS i_out = gain * i_ctrl

"""
@dataclass
class VCVS:
    name: str
    out_p: str
    out_n: str
    ctrl_p: str
    ctrl_n: str
    gain: float


@dataclass
class VCCS:
    name: str
    out_p: str
    out_n: str
    ctrl_p: str
    ctrl_n: str
    gain: float


@dataclass
class CCCS:
    name: str
    out_p: str
    out_n: str
    control_branch: str
    gain: float


@dataclass
class CCVS:
    name: str
    out_p: str
    out_n: str
    control_branch: str
    gain: float