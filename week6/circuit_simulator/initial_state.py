"""Jianhao Wu -- consistent initial states for small linear RC circuits.

Uses the existing Week 5 components and MNAAssembler without changing them.
Supported elements: Resistor, Capacitor, VoltageSource and CurrentSource,
with constant independent sources. Every capacitor needs an initial voltage,
measured from node_plus to node_minus. Inductors and controlled sources are
outside the scope of this worked example.

Construction temporarily introduces capacitor currents as extra unknowns:

    [G   B] [x0] = [b ]
    [B.T 0] [ic]   [vc]

B is the capacitor incidence matrix. This enforces KCL, source constraints
and the specified capacitor voltages together. The temporary currents are
not added to NodeMap; the returned x0 retains the existing MNA ordering.
The augmented system must be nonsingular (no redundant voltage constraints
or undetermined states). This is not a general DAE initialization solver.
"""

import numpy as np

from capacitor import Capacitor
from current_source import CurrentSource
from resistor import Resistor
from voltage_source import VoltageSource


def _rc_system(circuit, capacitor_voltages):
    """Assemble the existing circuit and validate the initialization inputs."""
    supported = (Resistor, Capacitor, VoltageSource, CurrentSource)
    for element in circuit.elements:
        if type(element) not in supported:
            raise ValueError(
                f"unsupported element {element.name!r}: "
                "initialization currently supports linear RC circuits "
                "with constant independent sources only"
            )

    capacitors = [e for e in circuit.elements if isinstance(e, Capacitor)]
    names = [e.name for e in capacitors]
    supplied = dict(capacitor_voltages)
    if set(supplied) != set(names):
        missing = sorted(set(names) - set(supplied))
        extra = sorted(set(supplied) - set(names))
        raise ValueError(
            f"provide exactly one voltage per capacitor; "
            f"missing={missing}, unknown={extra}"
        )
    vc = np.asarray([supplied[name] for name in names], dtype=float)
    if not np.all(np.isfinite(vc)):
        raise ValueError("initial capacitor voltages must be finite")

    G, M, b = circuit.assemble_system()
    G, M = G.toarray(), M.toarray()
    b = np.asarray(b, dtype=float)
    if b.size == 0:
        raise ValueError("the circuit has no MNA unknowns")
    if not all(np.all(np.isfinite(a)) for a in (G, M, b)):
        raise ValueError("assembled matrices and source vector must be finite")

    B = np.zeros((b.size, len(capacitors)))
    for column, capacitor in enumerate(capacitors):
        plus = circuit.node_map.get_index(capacitor.node_plus)
        minus = circuit.node_map.get_index(capacitor.node_minus)
        if plus is not None:
            B[plus, column] += 1.0
        if minus is not None:
            B[minus, column] -= 1.0

    # These rows enforce constant voltage sources: A_source @ xdot = 0.
    branch_rows = [
        circuit.node_map.get_branch_index(e.name)
        for e in circuit.elements if isinstance(e, VoltageSource)
    ]
    return G, M, b, B, vc, names, branch_rows


def build_initial_state(circuit, capacitor_voltages):
    """Return a complete x0 in the existing NodeMap ordering.

    Example: build_initial_state(circuit, {"C1": 2.0}).
    All elements must already have been added. The function assembles the
    circuit, so branch indices are available afterwards. Build x0 again if
    the circuit or its sources change. Negative capacitor voltages are valid.
    """
    G, M, b, B, vc, names, branch_rows = _rc_system(
        circuit, capacitor_voltages
    )
    n, count = B.shape
    augmented = np.block([[G, B], [B.T, np.zeros((count, count))]])
    rhs = np.concatenate((b, vc))

    # Row scaling helps balance KCL rows (amperes) and voltage constraints.
    scale = np.max(np.abs(augmented), axis=1)
    scale[scale == 0.0] = 1.0
    try:
        solution = np.linalg.solve(augmented / scale[:, None], rhs / scale)
    except np.linalg.LinAlgError as error:
        raise ValueError(
            "cannot determine a unique initial state: check floating nodes, "
            "conflicting or redundant capacitor/source voltage constraints; "
            "a general DAE initializer is required for unsupported topologies"
        ) from error

    x0 = solution[:n]
    report = check_initial_state(circuit, x0, capacitor_voltages)
    if not report["consistent"]:
        raise ValueError("initial-state check failed: " + "; ".join(report["issues"]))
    return x0


def check_initial_state(
    circuit, x0, capacitor_voltages, *,
    voltage_atol=1e-9, current_atol=1e-12,
    voltage_rate_atol=1e-9, rtol=1e-9,
):
    """Report whether x0 satisfies the RC initial-state constraints.

    Checks the specified capacitor voltages and whether a derivative exists
    such that M @ xdot + G @ x0 = b. Constant voltage sources also require
    A_source @ xdot = 0. A transient initial state need not satisfy G @ x0=b.

    Absolute tolerances are in volts, amperes and volts/second respectively.
    Invalid inputs raise ValueError; a physically inconsistent candidate
    returns a report with consistent=False and human-readable issues.
    """
    for value in (voltage_atol, current_atol, voltage_rate_atol, rtol):
        if not np.isfinite(value) or value < 0:
            raise ValueError("tolerances must be finite and nonnegative")

    G, M, b, B, vc, names, branch_rows = _rc_system(
        circuit, capacitor_voltages
    )
    x0 = np.asarray(x0, dtype=float)
    if x0.shape != b.shape or not np.all(np.isfinite(x0)):
        raise ValueError("x0 must be a finite one-dimensional MNA vector")

    actual_vc = B.T @ x0
    cap_errors = actual_vc - vc
    cap_limits = voltage_atol + rtol * np.maximum(np.abs(actual_vc), np.abs(vc))

    # Include differentiated constant-source constraints, so a checker cannot
    # conceal an incorrect source current by allowing its voltage to change.
    A_source = G[branch_rows, :]
    derivative_matrix = np.vstack((M, A_source))
    derivative_rhs = np.concatenate((b - G @ x0, np.zeros(len(branch_rows))))
    scale = np.max(np.abs(derivative_matrix), axis=1)
    scale[scale == 0.0] = 1.0
    xdot, _, _, _ = np.linalg.lstsq(
        derivative_matrix / scale[:, None], derivative_rhs / scale, rcond=None
    )

    residual = M @ xdot + G @ x0 - b
    row_atol = np.full(b.size, current_atol)
    row_atol[branch_rows] = voltage_atol
    row_limits = row_atol + rtol * (
        np.abs(M) @ np.abs(xdot) + np.abs(G) @ np.abs(x0) + np.abs(b)
    )
    source_rate_residual = A_source @ xdot
    source_rate_limits = voltage_rate_atol + rtol * (np.abs(A_source) @ np.abs(xdot))

    issues = []
    for name, error, limit in zip(names, cap_errors, cap_limits):
        if abs(error) > limit:
            issues.append(f"capacitor {name}: initial voltage error {error:.6e} V")
    for row in np.flatnonzero(np.abs(residual) > row_limits):
        unit = "V" if row in branch_rows else "A"
        issues.append(f"MNA row {row}: residual {residual[row]:.6e} {unit}")
    if np.any(np.abs(source_rate_residual) > source_rate_limits):
        issues.append("derivative violates a constant voltage-source constraint")
    if not np.all(np.isfinite(xdot)) or not np.all(np.isfinite(residual)):
        issues.append("initial-state calculation produced non-finite values")

    return {
        "consistent": not issues,
        "issues": issues,
        "capacitor_voltages": dict(zip(names, actual_vc.tolist())),
        "capacitor_voltage_errors": dict(zip(names, cap_errors.tolist())),
        "xdot": xdot,
        "mna_residual": residual,
        "source_voltage_rate_residual": source_rate_residual,
    }
