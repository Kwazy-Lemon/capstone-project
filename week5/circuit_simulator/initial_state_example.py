"""Jianhao Wu -- worked examples of complete, consistent RLC and RC initial states.

Place this file and initial_state.py in week5/circuit_simulator alongside
the existing components. Run: python initial_state_example.py

Dependencies: numpy, scipy, matplotlib.
Main example: a series inductor feeding a loaded two-stage capacitor ladder.
Example 1: a single capacitor starting at 2 V rather than its DC equilibrium.
Example 2: a loaded two-stage RC ladder with independently specified capacitor
voltages. It checks both node current balances, compares backward Euler with
an independent continuous-time solution and saves a figure plus waveform CSV
alongside this script. No existing solver or component is changed.
Parameter sweep: 2 circuits x 5 initial conditions x 5 time steps = 50 cases.
Run only the matrix with: python initial_state_example.py --sweep-only
The sweep saves per-case CSV results, a JSON summary and a convergence plot.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from capacitor import Capacitor
from inductor import Inductor
from initial_state import build_initial_state, check_initial_state
from linear_solvers import ac_sweep, backward_euler, dc_solve
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


def print_state(circuit, x0):
    """Display labels in the ordering actually assigned by the existing code."""
    print("Complete initial MNA state (ground is omitted):")
    for node, index in sorted(circuit.node_map.nodes.items(), key=lambda item: item[1]):
        print(f"  x0[{index}] = V({node}) = {x0[index]:.6f} V")
    for branch, index in sorted(circuit.node_map.branches.items(), key=lambda item: item[1]):
        print(f"  x0[{index}] = I({branch}) = {x0[index]:.6f} A")
    print("  I(Vs) is positive from the source's positive to negative terminal.")


def simple_rc_example():
    print("=== Example 1: single-capacitor RC circuit ===")
    source_voltage = 10.0
    resistance = 1000.0
    capacitance = 1e-6
    initial_voltage = 2.0

    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "in", "gnd", source_voltage))
    circuit.add_element(Resistor("R1", "in", "out", resistance))
    circuit.add_element(Capacitor("C1", "out", "gnd", capacitance))

    capacitor_voltages = {"C1": initial_voltage}
    x0 = build_initial_state(circuit, capacitor_voltages)
    report = check_initial_state(circuit, x0, capacitor_voltages)

    node_in = circuit.node_map.get_index("in")
    node_out = circuit.node_map.get_index("out")
    source_index = circuit.node_map.get_branch_index("Vs")

    print_state(circuit, x0)

    expected_current = (source_voltage - initial_voltage) / resistance
    expected = x0.copy()
    expected[node_in] = source_voltage
    expected[node_out] = initial_voltage
    expected[source_index] = -expected_current
    if not np.allclose(x0, expected, rtol=1e-10, atol=1e-12):
        raise RuntimeError("initial state does not match the independent RC calculation")
    if not report["consistent"]:
        raise RuntimeError("valid initial state was rejected")

    print(f"\nResistor current at t=0: {expected_current:.6f} A")
    print(f"Capacitor voltage rate at t=0: {report['xdot'][node_out]:.6f} V/s")
    print("Initial-state constraint check: PASS")

    # A supplied capacitor voltage alone is not a complete MNA state.
    # Setting the source current to zero would violate KCL at the input node.
    wrong_current = x0.copy()
    wrong_current[source_index] = 0.0
    bad_report = check_initial_state(circuit, wrong_current, capacitor_voltages)
    if bad_report["consistent"]:
        raise RuntimeError("incorrect source current was not detected")
    print("\nIncorrect source-current example: correctly rejected")
    for issue in bad_report["issues"]:
        print(f"  {issue}")

    # Check the stated capacitor voltage separately from the MNA equations.
    wrong_voltage = x0.copy()
    wrong_voltage[node_out] = 3.0
    wrong_voltage[source_index] = -(source_voltage - 3.0) / resistance
    voltage_report = check_initial_state(circuit, wrong_voltage, capacitor_voltages)
    if voltage_report["consistent"]:
        raise RuntimeError("incorrect capacitor initial voltage was not detected")
    print("\nIncorrect capacitor-voltage example: correctly rejected")
    for issue in voltage_report["issues"]:
        print(f"  {issue}")

    # Reuse the existing solver, rather than implementing transient analysis again.
    # This confirms the returned vector can be passed directly to backward_euler.
    G, M, b = circuit.assemble_system()
    tau = resistance * capacitance
    dt = tau / 100
    t, X = backward_euler(G, M, lambda time: b, x0, dt, 5 * tau)
    analytical = source_voltage + (initial_voltage - source_voltage) * np.exp(-t / tau)
    error = X[:, node_out] - analytical
    if not np.array_equal(X[0], x0):
        raise RuntimeError("transient solver did not preserve the supplied t=0 state")
    if not np.allclose(X[:, node_out], analytical, rtol=0.0, atol=0.02):
        raise RuntimeError("transient result disagrees with the RC analytical solution")

    print("\nExisting backward Euler solver: PASS")
    print(f"  dt = {dt:.6e} s, tstop = {t[-1]:.6e} s")
    print(f"  Maximum sampled voltage error: {np.max(np.abs(error)):.6e} V")
    print(f"  RMS sampled voltage error: {np.sqrt(np.mean(error ** 2)):.6e} V")
    print("Single-capacitor example: PASS")


def save_ladder_outputs(t, voltages, source_current, reference, source_reference):
    """Save artifacts suitable for displaying in Colab or at the group meeting."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    folder = Path(__file__).resolve().parent
    csv_path = folder / "initial_state_ladder_waveforms.csv"
    np.savetxt(
        csv_path,
        np.column_stack((t, voltages, source_current, reference, source_reference)),
        delimiter=",",
        header="time_s,VA_BE_V,VB_BE_V,IVs_BE_A,VA_exact_V,VB_exact_V,IVs_exact_A",
        comments="",
    )

    fig, axes = plt.subplots(
        3, 1, figsize=(11, 9), gridspec_kw={"height_ratios": [1, 2, 1.5]},
        layout="constrained",
    )
    fig.suptitle("Consistent initialization of a loaded two-stage RC ladder", fontsize=15)
    axes[0].axis("off")
    axes[0].text(
        0.05, 0.95,
        "10 V -- R1 (1 kohm) -- A -- R2 (2 kohm) -- B\n"
        "                      |                  |\n"
        "                  C1 (1 uF)        C2 (2 uF) || R3 (3 kohm)\n"
        "                      |                  |\n"
        "                     GND                GND\n"
        "Initial conditions: V(A) = 2 V, V(B) = 1 V, I(Vs) = -8 mA",
        transform=axes[0].transAxes, va="top", family="monospace", fontsize=10,
    )
    time_ms = t * 1000
    for column, label, color in ((0, "A", "#1768AC"), (1, "B", "#CB6015")):
        axes[1].plot(time_ms, voltages[:, column], color=color, label=f"V({label}), backward Euler")
        axes[1].plot(
            time_ms, reference[:, column], "--", color=color, alpha=0.7,
            label=f"V({label}), exact reference",
        )
        axes[1].scatter([0], [voltages[0, column]], color=color, zorder=3)
    axes[1].set(ylabel="Node voltage (V)", title="Two coupled capacitor voltages")
    axes[1].legend(loc="lower right", fontsize=9)
    axes[2].plot(time_ms, source_current * 1000, color="#1768AC", label="Backward Euler")
    axes[2].plot(time_ms, source_reference * 1000, "--", color="#CB6015", label="Exact reference")
    axes[2].set(xlabel="Time (ms)", ylabel="I(Vs) (mA)", title="Voltage-source branch current")
    axes[2].legend(fontsize=9)
    for ax in axes[1:]:
        ax.grid(alpha=0.25)
        ax.set_xlim(time_ms[0], time_ms[-1])
    image_path = folder / "initial_state_ladder.png"
    fig.savefig(image_path, dpi=160)
    plt.close(fig)
    print(f"Figure saved to {image_path}")
    print(f"Waveform data saved to {csv_path}")


def loaded_rc_example():
    print("\n=== Example 2: two-capacitor RC ladder with a resistive load ===")
    source_voltage = 10.0
    r1, r2, r3 = 1000.0, 2000.0, 3000.0
    c1, c2 = 1e-6, 2e-6
    va0, vb0 = 2.0, 1.0

    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "in", "gnd", source_voltage))
    circuit.add_element(Resistor("R1", "in", "A", r1))
    circuit.add_element(Resistor("R2", "A", "B", r2))
    circuit.add_element(Resistor("R3", "B", "gnd", r3))
    circuit.add_element(Capacitor("C1", "A", "gnd", c1))
    circuit.add_element(Capacitor("C2", "B", "gnd", c2))

    capacitor_voltages = {"C1": va0, "C2": vb0}
    x0 = build_initial_state(circuit, capacitor_voltages)
    report = check_initial_state(circuit, x0, capacitor_voltages)
    node_in = circuit.node_map.get_index("in")
    node_a = circuit.node_map.get_index("A")
    node_b = circuit.node_map.get_index("B")
    source_index = circuit.node_map.get_branch_index("Vs")
    print_state(circuit, x0)
    expected = np.empty_like(x0)
    expected[node_in], expected[node_a], expected[node_b] = source_voltage, va0, vb0
    expected[source_index] = -(source_voltage - va0) / r1
    if not report["consistent"] or not np.allclose(x0, expected, rtol=1e-10, atol=1e-12):
        raise RuntimeError("two-capacitor initial state failed its independent check")

    # Independent KCL calculation: all capacitor currents are toward ground.
    ir1 = (source_voltage - va0) / r1
    ir2 = (va0 - vb0) / r2
    ir3 = vb0 / r3
    ic1 = c1 * report["xdot"][node_a]
    ic2 = c2 * report["xdot"][node_b]
    kcl = np.array([ir1 - ir2 - ic1, ir2 - ir3 - ic2])
    expected_rates = np.array([(ir1 - ir2) / c1, (ir2 - ir3) / c2])
    if not np.allclose(kcl, 0, rtol=0, atol=1e-12):
        raise RuntimeError("node A or B violates initial KCL")
    if not np.allclose(report["xdot"][[node_a, node_b]], expected_rates, rtol=1e-10, atol=1e-8):
        raise RuntimeError("initial capacitor voltage rates are incorrect")
    print(f"\nResistor currents R1/R2/R3: {ir1:.9f}, {ir2:.9f}, {ir3:.9f} A")
    print(f"Capacitor currents C1/C2: {ic1:.9f}, {ic2:.9f} A")
    print(f"Initial dV(A)/dt: {expected_rates[0]:.6f} V/s")
    print(f"Initial dV(B)/dt: {expected_rates[1]:.6f} V/s")
    print(f"KCL residuals at A/B: {kcl[0]:.3e}, {kcl[1]:.3e} A")
    print("Two-capacitor initial-state and KCL checks: PASS")

    # The second capacitor adds an independent initial-voltage requirement.
    wrong_state = x0.copy()
    wrong_state[node_b] = 0.0
    wrong_report = check_initial_state(circuit, wrong_state, capacitor_voltages)
    if wrong_report["consistent"]:
        raise RuntimeError("incorrect C2 initial voltage was not detected")
    print("\nIncorrect C2 initial-voltage example: correctly rejected")
    for issue in wrong_report["issues"]:
        print(f"  {issue}")

    # Build a separate two-state reference directly from the circuit's KCL:
    # diag(C1,C2) vdot + K v = [Vs/R1, 0]. The matrix exponential gives the
    # continuous-time solution; it is independent of MNA time stepping.
    K = np.array([[1 / r1 + 1 / r2, -1 / r2], [-1 / r2, 1 / r2 + 1 / r3]])
    dynamics = -K / np.array([c1, c2])[:, None]
    # Independent DC divider values: 8.333333 V at A and 5 V at B.
    steady = source_voltage * np.array([r2 + r3, r3]) / (r1 + r2 + r3)
    initial = np.array([va0, vb0])
    G, M, b = circuit.assemble_system()
    dt, tstop = 1e-5, 0.025
    t, X = backward_euler(G, M, lambda time: b, x0, dt, tstop)
    reference = np.array([steady + expm(dynamics * time) @ (initial - steady) for time in t])
    source_reference = -(source_voltage - reference[:, 0]) / r1
    voltages = X[:, [node_a, node_b]]
    source_current = X[:, source_index]
    if not np.array_equal(X[0], x0):
        raise RuntimeError("transient solver changed the supplied initial state")
    if not np.allclose(voltages, reference, rtol=0, atol=0.02):
        raise RuntimeError("two-capacitor transient disagrees with its exact reference")
    if not np.allclose(source_current, source_reference, rtol=0, atol=2e-5):
        raise RuntimeError("source current disagrees with its exact reference")
    if not np.allclose(source_current + (source_voltage - voltages[:, 0]) / r1, 0, atol=1e-12, rtol=0):
        raise RuntimeError("input-node KCL failed during the transient")

    print("\nTwo-capacitor backward Euler / exact-reference comparison: PASS")
    print(f"  dt = {dt:.6e} s, tstop = {tstop:.6e} s")
    print(f"  DC equilibrium: V(A) = {steady[0]:.6f} V, V(B) = {steady[1]:.6f} V")
    errors = voltages - reference
    for column, label in enumerate(("A", "B")):
        maximum = np.max(np.abs(errors[:, column]))
        rms = np.sqrt(np.mean(errors[:, column] ** 2))
        print(f"  V({label}): maximum sampled error = {maximum:.6e} V, RMS = {rms:.6e} V")
    current_error = source_current - source_reference
    print(f"  I(Vs): maximum sampled error = {np.max(np.abs(current_error)):.6e} A")
    print(f"  I(Vs): RMS sampled error = {np.sqrt(np.mean(current_error ** 2)):.6e} A")
    save_ladder_outputs(t, voltages, source_current, reference, source_reference)


def sweep_circuit(topology, initial_voltages):
    """Build a sweep circuit, plus an independent reduced-state reference."""
    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vs", "in", "gnd", 10.0))
    if topology == "single_rc":
        circuit.add_element(Resistor("R1", "in", "out", 1000.0))
        circuit.add_element(Capacitor("C1", "out", "gnd", 1e-6))
        labels, capacitances = ["out"], np.array([1e-6])
        conductance = np.array([[1 / 1000.0]])
        steady = np.array([10.0])
    else:
        circuit.add_element(Resistor("R1", "in", "A", 1000.0))
        circuit.add_element(Resistor("R2", "A", "B", 2000.0))
        circuit.add_element(Resistor("R3", "B", "gnd", 3000.0))
        circuit.add_element(Capacitor("C1", "A", "gnd", 1e-6))
        circuit.add_element(Capacitor("C2", "B", "gnd", 2e-6))
        labels, capacitances = ["A", "B"], np.array([1e-6, 2e-6])
        conductance = np.array([[1 / 1000 + 1 / 2000, -1 / 2000],
                                [-1 / 2000, 1 / 2000 + 1 / 3000]])
        steady = np.array([10 * 5000 / 6000, 10 * 3000 / 6000])
    prescribed = {f"C{i + 1}": voltage for i, voltage in enumerate(initial_voltages)}
    return circuit, prescribed, labels, capacitances, conductance, steady


def save_sweep_plot(rows, folder):
    """Plot every initial-condition group; all signals use their own units."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), layout="constrained")
    for column, topology in enumerate(("single_rc", "loaded_rc")):
        for group in range(1, 6):
            selected = sorted(
                (r for r in rows if r["topology"] == topology and r["initial_group"] == group),
                key=lambda r: r["dt_us"],
            )
            if not all(np.isfinite(r.get("max_voltage_error_V", np.nan)) for r in selected):
                continue
            dt = [r["dt_us"] for r in selected]
            axes[0, column].loglog(dt, [r["max_voltage_error_V"] for r in selected], "o-", label=f"IC {group}")
            axes[1, column].loglog(dt, [r["rms_voltage_error_V"] for r in selected], "o-")
        axes[0, column].set_title("Single RC" if topology == "single_rc" else "Loaded two-stage RC")
        axes[0, column].set_ylabel("Worst-node maximum error (V)")
        axes[1, column].set(xlabel="Time step (microseconds)", ylabel="Worst-node RMS error (V)")
        axes[0, column].legend(fontsize=9)
        for ax in axes[:, column]:
            ax.grid(True, which="both", alpha=.25)
    fig.suptitle("50 distinct cases: 2 circuits x 5 initial conditions x 5 time steps")
    path = folder / "initial_state_sweep_convergence.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def run_parameter_sweep():
    """Run 50 distinct cases, with declared accuracy and convergence checks.

    These are deterministic parameter combinations, not repeated statistical
    trials. Component values and source voltage stay fixed in each topology.
    Accuracy limit: 2% of the initial voltage excursion, with a 1 V scale
    floor. KCL and voltage-source residual limits are 1e-10 A and 1e-9 V.
    A refinement passes only when both max and RMS errors, for every voltage
    and the source current, do not increase (relative slack 1e-6).
    """
    initial_sets = {
        "single_rc": [(0.0,), (2.0,), (5.0,), (-1.0,), (12.0,)],
        "loaded_rc": [(0.0, 0.0), (2.0, 1.0), (5.0, 3.0), (-1.0, 2.0), (12.0, 8.0)],
    }
    time_steps_us = [40.0, 20.0, 10.0, 5.0, 2.5]
    tstop = .025
    folder = Path(__file__).resolve().parent
    rows, reference_cache = [], {}
    print("\n=== Parameter sweep: 50 distinct circuit / initial-state / time-step cases ===")
    for topology, initial_cases in initial_sets.items():
        for group, initial in enumerate(initial_cases, start=1):
            previous_metrics = None
            for dt_us in time_steps_us:
                row = {"case_id": len(rows) + 1, "topology": topology,
                       "initial_group": group, "v1_initial_V": initial[0],
                       "v2_initial_V": initial[1] if len(initial) == 2 else "",
                       "dt_us": dt_us, "tstop_s": tstop, "passed": False,
                       "convergence": "BASELINE" if previous_metrics is None else "FAIL",
                       "max_voltage_error_V": np.nan, "rms_voltage_error_V": np.nan,
                       "max_source_error_A": np.nan, "rms_source_error_A": np.nan,
                       "max_kcl_residual_A": np.nan, "max_source_constraint_error_V": np.nan,
                       "voltage_error_limit_V": np.nan, "initial_state_pass": False,
                       "reference_accuracy_pass": False, "equations_pass": False,
                       "error": ""}
                try:
                    circuit, prescribed, labels, caps, K, steady = sweep_circuit(topology, initial)
                    x0 = build_initial_state(circuit, prescribed)
                    check = check_initial_state(circuit, x0, prescribed)
                    G, M, b = circuit.assemble_system()
                    nodes = [circuit.node_map.get_index(name) for name in labels]
                    source = circuit.node_map.get_branch_index("Vs")
                    input_node = circuit.node_map.get_index("in")
                    expected = x0.copy()
                    expected[input_node], expected[nodes] = 10.0, initial
                    expected[source] = -(10.0 - initial[0]) / 1000.0
                    t, X = backward_euler(G, M, lambda time: b, x0, dt_us * 1e-6, tstop)
                    row["initial_state_pass"] = bool(check["consistent"] and
                        np.allclose(x0, expected, rtol=0, atol=1e-10) and np.array_equal(X[0], x0))

                    key = (topology, dt_us)
                    if key not in reference_cache:
                        dynamics = -K / caps[:, None]
                        reference_cache[key] = np.array([expm(dynamics * time) for time in t])
                    reference = steady + np.einsum("tij,j->ti", reference_cache[key], np.array(initial) - steady)
                    source_reference = -(10.0 - reference[:, 0]) / 1000.0
                    voltage_error = X[:, nodes] - reference
                    source_error = X[:, source] - source_reference
                    maxima = np.max(np.abs(voltage_error), axis=0)
                    rms_values = np.sqrt(np.mean(voltage_error ** 2, axis=0))
                    metrics = np.concatenate((maxima, rms_values,
                        [np.max(np.abs(source_error)), np.sqrt(np.mean(source_error ** 2))]))
                    row["max_voltage_error_V"] = float(np.max(maxima))
                    row["rms_voltage_error_V"] = float(np.max(rms_values))
                    row["max_source_error_A"], row["rms_source_error_A"] = map(float, metrics[-2:])
                    limit = .02 * max(1.0, float(np.max(np.abs(np.array(initial) - steady))))
                    row["voltage_error_limit_V"] = limit
                    row["reference_accuracy_pass"] = bool(np.all(np.isfinite(X)) and
                        np.all(maxima <= limit) and metrics[-2] <= limit / 1000.0)

                    # Independent, directly evaluated node equations over every BE step.
                    v = X[:, nodes]
                    kcl = np.diff(v, axis=0) / (dt_us * 1e-6) * caps + v[1:] @ K.T
                    kcl[:, 0] -= 10.0 / 1000.0
                    source_kcl = X[:, source] + (X[:, input_node] - v[:, 0]) / 1000.0
                    row["max_kcl_residual_A"] = float(max(np.max(np.abs(kcl)), np.max(np.abs(source_kcl))))
                    row["max_source_constraint_error_V"] = float(np.max(np.abs(X[:, input_node] - 10.0)))
                    row["equations_pass"] = bool(row["max_kcl_residual_A"] <= 1e-10 and
                        row["max_source_constraint_error_V"] <= 1e-9)
                    refining = previous_metrics is None or np.all(metrics <= previous_metrics * (1 + 1e-6) + 1e-12)
                    row["convergence"] = "BASELINE" if previous_metrics is None else "PASS" if refining else "FAIL"
                    previous_metrics = metrics
                    row["passed"] = bool(row["initial_state_pass"] and row["reference_accuracy_pass"] and
                                         row["equations_pass"] and refining)
                except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                    row["error"] = str(error)
                    previous_metrics = None
                rows.append(row)
                print(f"{row['case_id']:02d}/50 {topology:9s} IC{group} dt={dt_us:4.1f} us "
                      f"max={row['max_voltage_error_V']:.6e} V  {'PASS' if row['passed'] else 'FAIL'}")

    csv_path = folder / "initial_state_sweep_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    groups = [{"topology": topology, "initial_group": group,
               "passed": all(r["passed"] for r in rows if r["topology"] == topology and r["initial_group"] == group)}
              for topology in initial_sets for group in range(1, 6)]
    summary = {"total_cases": len(rows), "passed_cases": sum(r["passed"] for r in rows),
               "failed_cases": sum(not r["passed"] for r in rows),
               "circuits": list(initial_sets), "initial_conditions": initial_sets,
               "time_steps_us": time_steps_us, "tstop_s": tstop,
               "criteria": {"reference_accuracy": "max error <= 2% of max(1 V, initial voltage excursion)",
                            "kcl_atol_A": 1e-10, "source_voltage_atol_V": 1e-9,
                            "convergence": "per-signal max and RMS errors do not increase when dt halves"},
               "refinement_groups": groups}
    json_path = folder / "initial_state_sweep_summary.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    plot_path = save_sweep_plot(rows, folder)
    print(f"\nPARAMETER SWEEP: {summary['passed_cases']}/{len(rows)} PASS; {summary['failed_cases']} FAIL")
    print(f"REFINEMENT GROUPS: {sum(g['passed'] for g in groups)}/{len(groups)} PASS")
    print("Coverage: 2 circuits x 5 initial conditions x 5 time steps; 25 ms per case")
    print("Time steps: 40, 20, 10, 5, 2.5 microseconds")
    print("Accuracy criterion: maximum voltage error <= 2% of max(1 V, initial voltage excursion)")
    print("Topology    dt(us)  worst max error(V)  worst RMS error(V)")
    for topology in initial_sets:
        for dt in time_steps_us:
            subset = [r for r in rows if r["topology"] == topology and r["dt_us"] == dt]
            print(f"{topology:9s} {dt:6.1f}  {max(r['max_voltage_error_V'] for r in subset):.6e}       "
                  f"{max(r['rms_voltage_error_V'] for r in subset):.6e}")
    print(f"Results: {csv_path.name}; {json_path.name}; {plot_path.name}")
    if summary["failed_cases"]:
        raise RuntimeError("parameter sweep contains failures; inspect initial_state_sweep_results.csv")
    return summary


def rlc_example():
    """Check a two-capacitor ladder with nonzero initial inductor current."""
    print("=== Main example: loaded two-stage RLC ladder ===")
    vs, r1, r2, r3 = 10.0, 100.0, 2000.0, 3000.0
    l1, c1, c2 = 0.04, 1e-6, 2e-6
    circuit = MNAAssembler()
    for element in (
        VoltageSource("Vs", "in", "0", vs),
        Resistor("R1", "in", "p", r1),
        Inductor("L1", "p", "A", l1),
        Capacitor("C1", "A", "0", c1),
        Resistor("R2", "A", "B", r2),
        Capacitor("C2", "B", "0", c2),
        Resistor("R3", "B", "0", r3),
    ):
        circuit.add_element(element)

    vc, il = {"C1": 2.0, "C2": 1.0}, {"L1": 0.004}
    x0 = build_initial_state(circuit, vc, inductor_currents=il)
    report = check_initial_state(circuit, x0, vc, inductor_currents=il)
    print_state(circuit, x0)
    np.testing.assert_allclose(x0, [10.0, 9.6, 2.0, 1.0, -0.004, 0.004])
    g, m, b = circuit.assemble_system()
    node = circuit.node_map.get_index
    branch = circuit.node_map.get_branch_index
    indices = [node("A"), node("B"), branch("L1")]
    np.testing.assert_allclose(report["xdot"][indices], [3500.0, 250.0 / 3.0, 190.0])
    wrong = x0.copy()
    wrong[branch("L1")] = 0.0
    if check_initial_state(circuit, wrong, vc, inductor_currents=il)["consistent"]:
        raise AssertionError("incorrect inductor initial current was accepted")

    # Independent ODE in [VA, VB, IL], derived from KCL and KVL.
    a = np.array([
        [-1 / (r2 * c1), 1 / (r2 * c1), 1 / c1],
        [1 / (r2 * c2), -(1 / r2 + 1 / r3) / c2, 0],
        [-1 / l1, 0, -r1 / l1],
    ])
    steady = np.linalg.solve(a, [0.0, 0.0, -vs / l1])
    dc = dc_solve(g, b)
    np.testing.assert_allclose(dc[indices], steady)
    previous_error = None
    for dt in (20e-6, 10e-6, 5e-6, 2.5e-6, 1.25e-6):
        t, x = backward_euler(g, m, lambda time: b, x0, dt, 0.025)
        transition = expm(a * dt)
        reference = np.empty((len(t), 3))
        reference[0] = x0[indices]
        for k in range(1, len(t)):
            reference[k] = steady + transition @ (reference[k - 1] - steady)
        error = x[:, indices] - reference
        metrics = np.stack((np.max(np.abs(error), axis=0), np.sqrt(np.mean(error**2, axis=0))))
        if previous_error is not None and not np.all(metrics < previous_error):
            raise AssertionError("RLC errors did not decrease under step refinement")
        previous_error = metrics
        print(f"  dt={dt * 1e6:g} us | max [V, V, A]: {metrics[0]} | RMS: {metrics[1]}")

    # Unit small-signal excitation and an independent complex impedance check.
    frequency = np.geomspace(10.0, 100_000.0, 401)
    excitation = np.zeros_like(b)
    excitation[branch("Vs")] = 1.0
    phasors = ac_sweep(g, m, excitation, frequency)
    s = 2j * np.pi * frequency
    z_b = 1 / (1 / r3 + s * c2)
    z_a = 1 / (s * c1 + 1 / (r2 + z_b))
    exact_gain = z_a / (r1 + s * l1 + z_a) * z_b / (r2 + z_b)
    gain = phasors[:, node("B")]
    np.testing.assert_allclose(gain, exact_gain, rtol=1e-10, atol=1e-12)
    print(f"  AC complex gain matches impedance reference at {len(frequency)} frequencies")
    print(f"  Peak |V(B)/V(in)| = {max(abs(gain)):.4f} at {frequency[np.argmax(abs(gain))]:.1f} Hz")
    print("RLC initialization, DC, AC, and transient checks passed.")

    import matplotlib.pyplot as plt

    folder = Path(__file__).resolve().parent
    for filename, columns, scale, ylabel, limit in (
        ("initial_state_rlc_voltages.png", (0, 1), 1, "Node voltage (V)", 6),
        ("initial_state_rlc_current.png", (2,), 1000, "Inductor current (mA)", 6),
    ):
        fig, ax = plt.subplots(figsize=(8, 3))
        for j in columns:
            label = ("V(A)", "V(B)", "I(L1)")[j]
            line, = ax.plot(t * 1000, x[:, indices[j]] * scale, label=f"{label}, Backward Euler")
            ax.plot(t * 1000, reference[:, j] * scale, "--", color=line.get_color(), alpha=0.7,
                    label=f"{label}, reference")
        ax.set(xlabel="Time (ms)", ylabel=ylabel, xlim=(0, limit))
        ax.grid(alpha=0.25)
        ax.legend()
        fig.savefig(folder / filename, bbox_inches="tight", dpi=160)
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.semilogx(frequency, 20 * np.log10(np.abs(gain)))
    ax.set(xlabel="Frequency (Hz)", ylabel="Output gain (dB)")
    ax.grid(which="both", alpha=0.25)
    fig.savefig(folder / "initial_state_rlc_ac.png", bbox_inches="tight", dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep-only", action="store_true", help="run the 50-case matrix without the worked demonstrations")
    parser.add_argument("--rlc-only", action="store_true", help="run only the RLC demonstration")
    args = parser.parse_args()
    if args.rlc_only and args.sweep_only:
        parser.error("choose either --rlc-only or --sweep-only")
    if not args.sweep_only:
        rlc_example()
    if args.rlc_only:
        return
    if not args.sweep_only:
        simple_rc_example()
        loaded_rc_example()
    run_parameter_sweep()
    print("\nAll initial-state example checks passed.")


if __name__ == "__main__":
    main()
