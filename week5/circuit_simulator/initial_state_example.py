"""Jianhao Wu -- worked examples of complete, consistent RC initial states.

Place this file and initial_state.py in week5/circuit_simulator alongside
the existing components. Run: python initial_state_example.py

Dependencies: numpy, scipy, matplotlib.
Example 1: a single capacitor starting at 2 V rather than its DC equilibrium.
Example 2: a loaded two-stage RC ladder with independently specified capacitor
voltages. It checks both node current balances, compares backward Euler with
an independent continuous-time solution and saves a figure plus waveform CSV
alongside this script. No existing solver or component is changed.
"""

from pathlib import Path

import numpy as np
from scipy.linalg import expm

from capacitor import Capacitor
from initial_state import build_initial_state, check_initial_state
from linear_solvers import backward_euler
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


def main():
    simple_rc_example()
    loaded_rc_example()
    print("\nAll initial-state example checks passed.")


if __name__ == "__main__":
    main()
