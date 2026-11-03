
"""DC analysis: build G and b, then solve G*x = b.

Public function:
    G, b, x, unknowns = mna(elements)

Input: the element objects returned by your netlist_parser.
Output: Python lists; unknowns[i] tells you what x[i] represents.

Gx + Cx' + f(x) = b
1. x' = 0; f(x) = 0; Gx = b 
2. find Gx and b from the netlist
3. x = G^-1 * b <- figure out solve for x 
"""
import math
from elements import (
    Resistor, Capacitor, Inductor, VoltageSource, CurrentSource,
    VCVS, VCCS, CCVS, CCCS,
)


def mna(elements):
    """Return G, b, x, unknowns for the DC system of equations.
    """
    elements = list(elements)
    if not elements:
        raise ValueError("No elements to analyse.")

    def node_name(name):
        name = str(name).lower()
        if name in ("0", "gnd", "ground"):
            return "0"
        else:
            return name

    def terminals(element):
        if isinstance(element, (Resistor, Capacitor, Inductor,
                                VoltageSource, CurrentSource)):
            return node_name(element.node1), node_name(element.node2)
        return node_name(element.out_p), node_name(element.out_n) # else it maybe VCVS, VCCS, CCVS, CCCS

    node_index = {}
    element_names = set()
    supported = (Resistor, Capacitor, Inductor, VoltageSource, CurrentSource,
                 VCVS, VCCS, CCVS, CCCS)
    for element in elements:
        if not isinstance(element, supported):
            raise ValueError(f"Unsupported element type: {type(element).__name__}")
        key = element.name.upper()
        if key in element_names:
            raise ValueError(f"Duplicate element name: {element.name}")
        element_names.add(key)
        value = element.gain if isinstance(element, (VCVS, VCCS, CCVS, CCCS)) else element.value
        if not math.isfinite(value):
            raise ValueError(f"{element.name}: value must be finite.")
        if isinstance(element, (Resistor, Capacitor, Inductor)) and value <= 0:
            raise ValueError(f"{element.name}: R/C/L value must be positive.")
        nodes = list(terminals(element))
        if isinstance(element, (VCVS, VCCS)):
            nodes += [node_name(element.ctrl_p), node_name(element.ctrl_n)]
        for node in nodes:
            if node != "0" and node not in node_index:
                node_index[node] = len(node_index)

    control_names = set()
    for element in elements:
        if isinstance(element, (CCVS, CCCS)):
            key = element.control_branch.upper()
            if key not in element_names:
                raise ValueError(f"{element.name}: unknown control branch {element.control_branch}")
            control_names.add(key)
    branch_index = {}
    unknowns = [f"V({node})" for node in node_index]
    for element in elements:
        if (isinstance(element, (VoltageSource, Inductor, VCVS, CCVS))
                or element.name.upper() in control_names):
            branch_index[element.name.upper()] = len(unknowns)
            unknowns.append(f"I({element.name})")

    size = len(unknowns)
    G = [[0.0 for _ in range(size)] for _ in range(size)]
    b = [0.0 for _ in range(size)]

    def add_G(row, col, value):
        if row is not None and col is not None:
            G[row][col] += value

    def add_b(row, value):
        if row is not None:
            b[row] += value

    for element in elements:
        first, second = terminals(element)
        p = node_index.get(first)
        n = node_index.get(second)
        k = branch_index.get(element.name.upper())
        if isinstance(element, (VCVS, VCCS)):
            cp = node_index.get(node_name(element.ctrl_p))
            cn = node_index.get(node_name(element.ctrl_n))

        if k is not None:
            add_G(p, k, 1.0)
            add_G(n, k, -1.0)

            if isinstance(element, (VoltageSource, Inductor, VCVS, CCVS, Resistor)):
                add_G(k, p, 1.0)
                add_G(k, n, -1.0)
                if isinstance(element, VoltageSource):
                    add_b(k, element.value)            # vp - vn = V
                elif isinstance(element, Inductor):
                    pass                               # vp - vn = 0 in DC
                elif isinstance(element, VCVS):
                    add_G(k, cp, -element.gain)         # vp-vn-mu*(vcp-vcn)=0
                    add_G(k, cn, element.gain)
                elif isinstance(element, CCVS):
                    ctrl = branch_index[element.control_branch.upper()]
                    add_G(k, ctrl, -element.gain)       # vp-vn-r*i_control=0
                elif isinstance(element, Resistor):
                    add_G(k, k, -element.value)         # vp-vn-R*i=0
            else:
                add_G(k, k, 1.0)
                if isinstance(element, CurrentSource):
                    add_b(k, element.value)            # i = I
                elif isinstance(element, Capacitor):
                    pass                               # i = 0 in DC
                elif isinstance(element, VCCS):
                    add_G(k, cp, -element.gain)         # i-gm*(vcp-vcn)=0
                    add_G(k, cn, element.gain)
                elif isinstance(element, CCCS):
                    ctrl = branch_index[element.control_branch.upper()]
                    add_G(k, ctrl, -element.gain)       # i-beta*i_control=0

        elif isinstance(element, Resistor):
            conductance = 1.0 / element.value
            add_G(p, p, conductance)
            add_G(p, n, -conductance)
            add_G(n, p, -conductance)
            add_G(n, n, conductance)

        elif isinstance(element, CurrentSource):
            add_b(p, -element.value)
            add_b(n, element.value)

        elif isinstance(element, Capacitor):
            pass                                       # open circuit in DC

        elif isinstance(element, VCCS):
            add_G(p, cp, element.gain)
            add_G(p, cn, -element.gain)
            add_G(n, cp, -element.gain)
            add_G(n, cn, element.gain)

        elif isinstance(element, CCCS):
            ctrl = branch_index[element.control_branch.upper()]
            add_G(p, ctrl, element.gain)
            add_G(n, ctrl, -element.gain)

    x = _solve_linear(G, b)
    return G, b, x, unknowns


def _solve_linear(G, b):
    """
    G and b are preserved. Suitable for small coursework circuits, not large
    sparse or severely ill-conditioned systems. main.py only calls mna().
    """
    size = len(b)
    if size == 0:
        return []
    error = ("Singular or ill-conditioned DC system: no unique, numerically reliable solution. Check ground, "
             "floating nodes and conflicting/redundant ideal voltage constraints.")
    # Work on an augmented copy [G | b]; scale each equation first.
    a = []
    for row, rhs in zip(G, b):
        scale = max(abs(value) for value in row)
        if scale == 0.0:
            raise ValueError(error)
        a.append([value / scale for value in row] + [rhs / scale])

    # Forward elimination，convert G to an upper triangular matrix.
    for col in range(size):
        pivot = max(range(col, size), key=lambda row: abs(a[row][col]))
        """
        if abs(a[pivot][col]) < 1e-14:
            raise ValueError(error)
        """
        a[col], a[pivot] = a[pivot], a[col]
        for row in range(col + 1, size):
            factor = a[row][col] / a[col][col]
            a[row][col] = 0.0
            for j in range(col + 1, size + 1):
                a[row][j] -= factor * a[col][j]

    # Back substitution，compute the last unknown first, then work upwards.
    x = [0.0 for _ in range(size)]
    for row in range(size - 1, -1, -1):
        known = sum(a[row][j] * x[j] for j in range(row + 1, size))
        x[row] = (a[row][size] - known) / a[row][row]
    if not all(math.isfinite(value) for value in x):
        raise ValueError(error)
    return x
