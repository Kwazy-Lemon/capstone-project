"""Build actual components, solve them, and export Week 3 validation evidence.

Run from the repository root: python week3/circuit_simulator/linear_examples.py
All matrices come from MNAAssembler; no circuit matrix is entered by hand.
The numerical solvers are the unchanged Week 2 implementation.
"""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Save figures without needing an interactive GUI window.
import matplotlib.pyplot as plt
import numpy as np

from controlled_sources import VCVS, VCCS, CCCS, CCVS
from inductor import Inductor
from linear_solvers import dc_solve, ac_sweep, backward_euler
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


def save_csv(path, fields, rows):
    """Keep column labels and units alongside the numerical evidence."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(rows)


def controlled_source_examples(output_dir):
    """Four separate circuits, all with a 2 V control and a 2 kohm load.

    A 1 kohm resistor and a series 0 V probe produce Icontrol = 2 mA.
    The probe does not change the control circuit's voltage, but it makes
    that current available as an MNA unknown for CCCS and CCVS.
    """
    cases = [
        (VCVS("E1", "out", "GND", "control", "GND", 3.0), 6.0),
        (VCCS("G1", "out", "GND", "control", "GND", 0.001), -4.0),
        (CCCS("F1", "out", "GND", "Vsense", 2.0), -8.0),
        (CCVS("H1", "out", "GND", "Vsense", 1000.0), 2.0),
    ]
    rows = []
    for source, expected_voltage in cases:
        circuit = MNAAssembler()
        # Add the controlled source first to demonstrate forward references:
        # Vsense is named before it has been registered with the assembler.
        circuit.add_element(source)
        circuit.add_element(VoltageSource("Vin", "control", "GND", 2.0))
        circuit.add_element(Resistor("Rcontrol", "control", "sense", 1000.0))
        circuit.add_element(VoltageSource("Vsense", "sense", "GND", 0.0))
        circuit.add_element(Resistor("Rload", "out", "GND", 2000.0))

        G, M, b = circuit.assemble_system()
        x = dc_solve(G, b)
        vout = x[circuit.node_map.get_index("out")]
        icontrol = x[circuit.node_map.get_branch_index("Vsense")]
        np.testing.assert_allclose(vout, expected_voltage, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(icontrol, 0.002, rtol=1e-12, atol=1e-14)
        np.testing.assert_allclose(G @ x - b, 0.0, atol=1e-12)
        assert M.nnz == 0  # These four ideal sources store no energy.

        # Negative VCCS/CCCS voltages are intentional. Their positive current
        # leaves out toward ground, so load KCL gives Vout/Rload + Iout = 0.
        kind = type(source).__name__
        error = abs(vout - expected_voltage)
        print(f"{kind}: Vout = {vout: .6f} V, expected = {expected_voltage: .6f} V")
        rows.append((kind, vout, expected_voltage, error, icontrol))

    save_csv(output_dir / "controlled_sources.csv",
             ["source_type", "computed_V", "expected_V", "abs_error_V", "control_current_A"],
             rows)
    return {"controlled_source_cases": len(rows), "controlled_source_checks": "passed"}


def rl_example(output_dir):
    """Series R-L: Vin -- R -- out -- L -- GND, driven by 5 V.

    The voltage across L is Vout; I_L points from out to ground.
    R = 100 ohm, L = 10 mH, so tau = L/R = 100 us and I_final = 50 mA.
    """
    resistance, inductance, voltage = 100.0, 0.01, 5.0
    tau = inductance / resistance
    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vin", "in", "GND", voltage))
    circuit.add_element(Resistor("R1", "in", "out", resistance))
    circuit.add_element(Inductor("L1", "out", "GND", inductance))
    G, M, b = circuit.assemble_system()

    # Use name lookups instead of assuming that a voltage/current has a fixed
    # index. The assembler may assign different indices if element order changes.
    node_in = circuit.node_map.get_index("in")
    node_out = circuit.node_map.get_index("out")
    current_l = circuit.node_map.get_branch_index("L1")
    current_v = circuit.node_map.get_branch_index("Vin")

    # DC: dI_L/dt = 0, so the inductor has zero terminal voltage.
    xdc = dc_solve(G, b)
    np.testing.assert_allclose(xdc[node_out], 0.0, atol=1e-12)
    np.testing.assert_allclose(xdc[current_l], voltage / resistance, rtol=1e-12)
    print(f"RL DC: V_L = {xdc[node_out]:.6f} V, I_L = {xdc[current_l]:.6f} A")

    # AC: use a 5 V real phasor as the excitation for this example.
    # This is an explicit AC amplitude choice, independent of a DC operating point.
    frequencies = np.logspace(1, 5, 101)
    ac_states = ac_sweep(G, M, b.astype(complex), frequencies)
    impedance = 1j * 2 * np.pi * frequencies * inductance
    expected_i = voltage / (resistance + impedance)
    expected_v = impedance * expected_i
    np.testing.assert_allclose(ac_states[:, current_l], expected_i, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(ac_states[:, node_out], expected_v, rtol=1e-12, atol=1e-12)
    ac_error = float(np.max(np.abs(ac_states[:, current_l] - expected_i)))
    save_csv(output_dir / "rl_ac.csv",
             ["frequency_Hz", "I_L_real_A", "I_L_imag_A", "abs_current_error_A"],
             zip(frequencies, ac_states[:, current_l].real,
                 ac_states[:, current_l].imag, np.abs(ac_states[:, current_l] - expected_i)))

    # At t=0+ the source has stepped to 5 V, but I_L is still zero. Therefore
    # the resistor drop and source current are zero and Vout is initially 5 V.
    # Setting every entry of x0 to zero would violate the source constraint.
    x0 = np.zeros(b.size)
    x0[node_in] = voltage
    x0[node_out] = voltage
    errors = []
    convergence_rows = []
    for steps_per_tau in (20, 40):
        dt = tau / steps_per_tau
        time, states = backward_euler(G, M, lambda t: b, x0, dt, 5 * tau)
        exact_i = voltage / resistance * (1 - np.exp(-time / tau))

        # Independently check the backward Euler recurrence. This separates
        # a component/sign error from the expected time-discretization error.
        step_numbers = np.arange(time.size)
        discrete_i = voltage / resistance * (1 - (1 + dt / tau) ** (-step_numbers))
        np.testing.assert_allclose(states[:, current_l], discrete_i, rtol=1e-11, atol=1e-13)
        np.testing.assert_allclose(states[:, current_v], -states[:, current_l], atol=1e-13)
        max_error = float(np.max(np.abs(states[:, current_l] - exact_i)))
        errors.append(max_error)
        convergence_rows.append((dt, time.size - 1, max_error))

    if errors[1] >= errors[0]:
        raise AssertionError("Halving the time step should reduce the RL current error.")
    save_csv(output_dir / "rl_timestep_error.csv", ["dt_s", "steps", "max_current_error_A"],
             convergence_rows)
    # The arrays below are from the finer of the two simulations.
    save_csv(output_dir / "rl_transient.csv",
             ["time_s", "I_L_A", "analytic_I_L_A", "V_L_V", "analytic_V_L_V"],
             zip(time, states[:, current_l], exact_i,
                 states[:, node_out], voltage * np.exp(-time / tau)))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(time * 1e6, exact_i * 1e3, label="Analytical", linewidth=2)
    axes[0].plot(time * 1e6, states[:, current_l] * 1e3, "--", label="Backward Euler")
    axes[0].set(xlabel="Time (us)", ylabel="Inductor current (mA)", title="RL step response")
    axes[1].semilogx(frequencies, np.abs(expected_v / voltage), label="Analytical", linewidth=2)
    axes[1].semilogx(frequencies[::5], np.abs(ac_states[::5, node_out] / voltage),
                     "o", markersize=3, label="Assembled MNA")
    axes[1].set(xlabel="Frequency (Hz)", ylabel="|V_L / V_in|", title="RL frequency response")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "rl_validation.png", dpi=160)
    plt.close(fig)

    print(f"RL AC: maximum current error = {ac_error:.3e} A")
    print(f"RL transient: halving dt reduces maximum error from {errors[0]:.3e} to {errors[1]:.3e} A")
    return {
        "rl_dc_check": "passed", "rl_ac_check": "passed",
        "rl_ac_max_current_error_A": ac_error,
        "rl_transient_discrete_check": "passed", "rl_timestep_convergence": "passed",
        "rl_transient_fine_max_error_A": errors[1],
        "rl_time_constant_s": tau,
    }


def main():
    output_dir = Path(__file__).resolve().parent / "results"
    output_dir.mkdir(exist_ok=True)
    metrics = controlled_source_examples(output_dir)
    metrics.update(rl_example(output_dir))
    metrics["reference_method"] = "Analytical formulas and backward Euler recurrence; MATLAB not run."
    (output_dir / "validation_summary.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    print(f"Validation evidence saved in {output_dir}")


if __name__ == "__main__":
    main()
