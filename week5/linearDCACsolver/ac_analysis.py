"""Single-frequency linear AC analysis; shares elements.py with the DC solver.

    A, b, x, names = solve_ac(elements, f_hz, {"Vin": 1 + 0j})

The phasor convention is exp(+j*omega*t). Voltages are relative to ground;
branch current is positive from the first terminal to the second. A and b
are returned before elimination. No file I/O and no third-party packages.
"""
import cmath
import math

from elements import (Resistor, Capacitor, Inductor, VoltageSource, CurrentSource,
                      VCVS, VCCS, CCVS, CCCS)


def _node(name):
    name = str(name).lower()
    return "0" if name == "gnd" else name


def _terminals(e):
    if isinstance(e, (VCVS, VCCS, CCVS, CCCS)):
        return _node(e.out_p), _node(e.out_n)
    return _node(e.node1), _node(e.node2)


def _indices(elements):
    """Number node voltages first, then branch currents, just as in DC."""
    node, by_name, sensed = {}, {}, set()
    for e in elements:
        if not isinstance(e, (Resistor, Capacitor, Inductor, VoltageSource,
                              CurrentSource, VCVS, VCCS, CCVS, CCCS)):
            raise ValueError(f"Unsupported element: {type(e).__name__}")
        key = e.name.upper()
        if key in by_name:
            raise ValueError(f"Duplicate element name: {e.name}")
        by_name[key] = e
        if isinstance(e, (Resistor, Capacitor, Inductor)):
            if not math.isfinite(e.value) or e.value <= 0:
                raise ValueError(f"{e.name}: R, L and C must be finite and positive.")
        elif isinstance(e, (VCVS, VCCS, CCVS, CCCS)) and not math.isfinite(e.gain):
            raise ValueError(f"{e.name}: gain must be finite.")
        terminals = _terminals(e)
        if isinstance(e, (VCVS, VCCS)):
            terminals += (_node(e.ctrl_p), _node(e.ctrl_n))
        for label in terminals:
            if label != "0" and label not in node:
                node[label] = len(node)
        if isinstance(e, (CCVS, CCCS)):
            sensed.add(e.control_branch.upper())

    missing = sensed.difference(by_name)
    if missing:
        raise ValueError("Unknown controlling branch: " + ", ".join(sorted(missing)))
    names = [f"V({label})" for label in node]
    branch = {}
    for key, e in by_name.items():
        if isinstance(e, (VoltageSource, Inductor, VCVS, CCVS)) or key in sensed:
            branch[key] = len(names)
            names.append(f"I({e.name})")
    return node, branch, names, by_name


def solve_ac(elements, f_hz, ac_sources):
    """Return (A, b, x, names) at f_hz > 0.

    ac_sources maps independent source names (case-insensitive) to complex
    voltages or currents. For example, 1j means magnitude 1 at +90 degrees.
    A source not listed has zero AC excitation, regardless of its DC value.
    """
    elements = list(elements)
    if not elements:
        raise ValueError("No elements to analyse.")
    if not math.isfinite(f_hz) or f_hz <= 0:
        raise ValueError("AC frequency must be finite and positive, in Hz.")
    jw = 1j * (math.tau * f_hz)
    node, branch, names, by_name = _indices(elements)
    source = {}
    for name, phasor in ac_sources.items():
        key = name.upper()
        if not isinstance(by_name.get(key), (VoltageSource, CurrentSource)):
            raise ValueError(f"{name}: AC excitation must name an independent source.")
        value = complex(phasor)
        if key in source or not cmath.isfinite(value):
            raise ValueError(f"{name}: duplicate or non-finite AC excitation.")
        source[key] = value

    size = len(names)
    A = [[0j] * size for _ in range(size)]
    b = [0j] * size

    def add(row, col, value):
        if row is not None and col is not None:  # ground has no row or column
            A[row][col] += value

    def rhs(row, value):
        if row is not None:
            b[row] += value

    def voltage(row, pos, neg, coeff=1):
        """Add coeff*(Vpos - Vneg) to this equation."""
        add(row, pos, coeff)
        add(row, neg, -coeff)

    for e in elements:
        pos, neg = _terminals(e)
        p, n = node.get(pos), node.get(neg)
        key = e.name.upper()
        k = branch.get(key)
        if isinstance(e, (VCVS, VCCS)):
            cp, cn = node.get(_node(e.ctrl_p)), node.get(_node(e.ctrl_n))
        if isinstance(e, (CCVS, CCCS)):
            ctrl = branch[e.control_branch.upper()]
        if k is not None:
            # The two node equations count this branch current as leaving p.
            add(p, k, 1)
            add(n, k, -1)

        if isinstance(e, (Resistor, Capacitor)):
            Y = 1 / e.value if isinstance(e, Resistor) else jw * e.value
            if k is None:  # node KCL only; current = Y*(Vp - Vn)
                voltage(p, p, n, Y)
                voltage(n, p, n, -Y)
            elif isinstance(e, Resistor):  # retained current: Vp - Vn - R*I = 0
                voltage(k, p, n)
                add(k, k, -e.value)
            else:  # retained capacitor current: I - jw*C*(Vp - Vn) = 0
                add(k, k, 1)
                voltage(k, p, n, -Y)
        elif isinstance(e, Inductor):
            voltage(k, p, n)
            add(k, k, -jw * e.value)  # Vp - Vn - jw*L*I = 0
        elif isinstance(e, VoltageSource):
            voltage(k, p, n)
            rhs(k, source.get(key, 0j))
        elif isinstance(e, CurrentSource):
            I = source.get(key, 0j)
            if k is None:
                rhs(p, -I)
                rhs(n, I)
            else:
                add(k, k, 1)
                rhs(k, I)
        elif isinstance(e, VCVS):
            voltage(k, p, n)
            voltage(k, cp, cn, -e.gain)
        elif isinstance(e, CCVS):
            voltage(k, p, n)
            add(k, ctrl, -e.gain)
        elif isinstance(e, VCCS):
            if k is None:
                voltage(p, cp, cn, e.gain)
                voltage(n, cp, cn, -e.gain)
            else:
                add(k, k, 1)
                voltage(k, cp, cn, -e.gain)
        elif isinstance(e, CCCS):
            if k is None:
                add(p, ctrl, e.gain)
                add(n, ctrl, -e.gain)
            else:
                add(k, k, 1)
                add(k, ctrl, -e.gain)

    return A, b, _solve_complex(A, b), names


def _solve_complex(A, b):
    """Gaussian elimination on a scaled copy; pivots use complex magnitude."""
    error = "Singular or ill-conditioned AC system; check ground and ideal sources."
    M = []
    for row, rhs in zip(A, b):
        scale = max(map(abs, row))
        if scale == 0 or not math.isfinite(scale):
            raise ValueError(error)
        M.append([v / scale for v in row] + [rhs / scale])
    size = len(b)
    for col in range(size):
        pivot = max(range(col, size), key=lambda r: abs(M[r][col]))
        if abs(M[pivot][col]) < 1e-14:
            raise ValueError(error)
        M[col], M[pivot] = M[pivot], M[col]
        for r in range(col + 1, size):
            factor = M[r][col] / M[col][col]
            if factor:
                for j in range(col + 1, size + 1):
                    M[r][j] -= factor * M[col][j]
            M[r][col] = 0j
    x = [0j] * size
    for r in range(size - 1, -1, -1):
        known = sum(M[r][j] * x[j] for j in range(r + 1, size))
        x[r] = (M[r][-1] - known) / M[r][r]
    if not all(cmath.isfinite(v) for v in x):
        raise ValueError(error)
    return x
