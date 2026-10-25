"""Validate the Week 4 capacitor and independent current source.

Run from the repository root: python week4/circuit_simulator/linear_examples.py
Each circuit is built from components; no MNA circuit matrix is entered by hand.
The sparse DC/AC/backward Euler solvers are reused unchanged from Week 3.
"""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Save figures directly; no desktop plotting window is needed.
import matplotlib.pyplot as plt
import numpy as np

from capacitor import Capacitor
from current_source import CurrentSource
from linear_solvers import dc_solve, ac_sweep, backward_euler
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


def save_csv(path, fields, rows):
    """Save numerical evidence with descriptive column names and units."""
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(rows)


def current_source_example(output_dir):
    """Inject current from ground into a node with a 1 kohm resistor to ground."""
    resistance = 1000.0
    rows = []
    for current in (0.001, 0.0, -0.001):
        circuit = MNAAssembler()
        # Current is directed GND -> out, so positive values inject into out.
        circuit.add_element(CurrentSource("I1", "GND", "out", current))
        circuit.add_element(Resistor("R1", "out", "GND", resistance))
        G, M, b = circuit.assemble_system()
        state = dc_solve(G, b)
        vout = state[circuit.node_map.get_index("out")]
        expected = current * resistance
        np.testing.assert_allclose(vout, expected, atol=1e-12)
        np.testing.assert_allclose(G @ state - b, 0.0, atol=1e-14)
        assert M.nnz == 0
        rows.append((current, vout, expected, abs(vout - expected)))
        print(f"Current source: I = {current:+.6f} A, Vout = {vout:+.6f} V")

    save_csv(output_dir / "current_source_dc.csv",
             ["source_current_A", "computed_V", "expected_V", "abs_error_V"], rows)
    return {"current_source_dc_cases": len(rows), "current_source_dc_check": "passed"}


def rc_example(output_dir):
    """Voltage-driven RC low-pass: Vin -- R -- out, with C from out to ground."""
    resistance, capacitance, voltage = 1000.0, 1e-6, 5.0
    tau = resistance * capacitance
    cutoff = 1 / (2 * np.pi * tau)
    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vin", "in", "GND", voltage))
    circuit.add_element(Resistor("R1", "in", "out", resistance))
    circuit.add_element(Capacitor("C1", "out", "GND", capacitance))
    G, M, b = circuit.assemble_system()
    node_in = circuit.node_map.get_index("in")
    node_out = circuit.node_map.get_index("out")
    current_v = circuit.node_map.get_branch_index("Vin")

    # At DC a capacitor draws zero current. With no resistor drop, Vout = Vin.
    dc_state = dc_solve(G, b)
    np.testing.assert_allclose(dc_state[node_out], voltage, atol=1e-12)
    np.testing.assert_allclose(dc_state[current_v], 0.0, atol=1e-14)
    print(f"RC DC: Vout = {dc_state[node_out]:.6f} V, I_source = {dc_state[current_v]:.6f} A")

    # The AC source amplitude is explicitly chosen as 1 V. The stored 5 V DC
    # value does not automatically define an AC phasor in a general simulator.
    frequencies = np.unique(np.r_[np.logspace(0, 5, 101), cutoff])
    b_ac = np.zeros(b.size, dtype=complex)
    b_ac[current_v] = 1.0
    ac_states = ac_sweep(G, M, b_ac, frequencies)
    gain = ac_states[:, node_out] / ac_states[:, node_in]
    expected_gain = 1 / (1 + 1j * 2 * np.pi * frequencies * tau)
    np.testing.assert_allclose(gain, expected_gain, rtol=1e-12, atol=1e-14)
    cutoff_index = np.argmin(np.abs(frequencies - cutoff))
    np.testing.assert_allclose(abs(gain[cutoff_index]), 1 / np.sqrt(2), rtol=1e-12)
    ac_error = float(np.max(np.abs(gain - expected_gain)))
    save_csv(output_dir / "rc_ac.csv",
             ["frequency_Hz", "gain_real", "gain_imag", "gain_dB", "phase_deg", "abs_error"],
             zip(frequencies, gain.real, gain.imag, 20 * np.log10(np.abs(gain)),
                 np.angle(gain, deg=True), np.abs(gain - expected_gain)))

    # Immediately after a 5 V step the capacitor remains at its initial 0 V.
    # Vin is already 5 V, and the ideal source supplies 5/R amperes. Its MNA
    # branch current is negative because its reference direction is into ground.
    x0 = np.zeros(b.size)
    x0[node_in] = voltage
    x0[current_v] = -voltage / resistance
    convergence_rows = []
    for steps_per_tau in (20, 40):
        dt = tau / steps_per_tau
        time, states = backward_euler(G, M, lambda t: b, x0, dt, 5 * tau)
        exact_voltage = voltage * (1 - np.exp(-time / tau))

        # Deriving backward Euler directly from RC*dv/dt + v = Vin gives
        # v[n] = Vin*(1-q**n), q = 1/(1+dt/tau). This is an independent
        # discrete reference; the continuous curve still has timestep error.
        q = 1 / (1 + dt / tau)
        discrete_voltage = voltage * (1 - q ** np.arange(time.size))
        np.testing.assert_allclose(states[:, node_out], discrete_voltage, atol=1e-12)
        resistor_current = (voltage - states[1:, node_out]) / resistance
        capacitor_current = capacitance * np.diff(states[:, node_out]) / dt
        np.testing.assert_allclose(capacitor_current, resistor_current, atol=1e-13)
        max_error = float(np.max(np.abs(states[:, node_out] - exact_voltage)))
        convergence_rows.append((dt, time.size - 1, max_error))

    coarse_error, fine_error = convergence_rows[0][2], convergence_rows[1][2]
    if fine_error >= coarse_error:
        raise AssertionError("Halving the RC time step should reduce voltage error.")
    # These arrays retain the finer simulation from the loop above.
    save_csv(output_dir / "rc_transient.csv",
             ["time_s", "computed_voltage_V", "analytic_voltage_V", "abs_error_V"],
             zip(time, states[:, node_out], exact_voltage,
                 np.abs(states[:, node_out] - exact_voltage)))
    save_csv(output_dir / "rc_timestep_error.csv", ["dt_s", "steps", "max_voltage_error_V"],
             convergence_rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(time * 1e3, exact_voltage, label="Analytical", linewidth=2)
    axes[0].plot(time * 1e3, states[:, node_out], "--", label="Backward Euler")
    axes[0].set(xlabel="Time (ms)", ylabel="Capacitor voltage (V)", title="RC step response")
    axes[1].semilogx(frequencies, 20 * np.log10(np.abs(expected_gain)),
                     label="Analytical", linewidth=2)
    axes[1].semilogx(frequencies[::5], 20 * np.log10(np.abs(gain[::5])),
                     "o", markersize=3, label="Assembled MNA")
    axes[1].set(xlabel="Frequency (Hz)", ylabel="Gain (dB)", title="RC low-pass response")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "rc_validation.png", dpi=160)
    plt.close(fig)

    print(f"RC AC: cutoff = {cutoff:.6f} Hz, maximum gain error = {ac_error:.3e}")
    print(f"RC transient: halving dt reduces error from {coarse_error:.3e} to {fine_error:.3e} V")
    return {
        "capacitor_dc_check": "passed", "rc_ac_check": "passed",
        "rc_ac_max_gain_error": ac_error, "rc_cutoff_Hz": cutoff,
        "rc_transient_discrete_check": "passed", "rc_timestep_convergence": "passed",
        "rc_transient_coarse_max_error_V": coarse_error,
        "rc_transient_fine_max_error_V": fine_error,
    }


def constant_current_charge_example(output_dir):
    """A 1 mA source charges a 1 uF capacitor starting at 0.25 V.

    C*dv/dt = I gives v(t) = v0 + I*t/C. There is no DC equilibrium for this
    circuit with nonzero I; only its well-defined transient is solved here.
    """
    current, capacitance, initial_voltage = 0.001, 1e-6, 0.25
    circuit = MNAAssembler()
    circuit.add_element(CurrentSource("I1", "GND", "out", current))
    circuit.add_element(Capacitor("C1", "out", "GND", capacitance))
    G, M, b = circuit.assemble_system()
    node_out = circuit.node_map.get_index("out")
    x0 = np.zeros(b.size)
    x0[node_out] = initial_voltage
    time, states = backward_euler(G, M, lambda t: b, x0, 10e-6, 1e-3)
    expected = initial_voltage + current * time / capacitance
    actual = states[:, node_out]
    # A linear ramp is integrated exactly by backward Euler up to roundoff.
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    max_error = float(np.max(np.abs(actual - expected)))
    save_csv(output_dir / "constant_current_charge.csv",
             ["time_s", "computed_voltage_V", "analytic_voltage_V", "abs_error_V"],
             zip(time, actual, expected, np.abs(actual - expected)))

    fig, axis = plt.subplots(figsize=(6, 4))
    axis.plot(time * 1e3, expected, label="Analytical: v0 + I*t/C", linewidth=2)
    axis.plot(time[::5] * 1e3, actual[::5], "o", markersize=3, label="Backward Euler")
    axis.set(xlabel="Time (ms)", ylabel="Capacitor voltage (V)", title="Constant-current charging")
    axis.grid(alpha=0.25)
    axis.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "constant_current_charge.png", dpi=160)
    plt.close(fig)

    print(f"Constant-current charge: final voltage = {actual[-1]:.6f} V, error = {max_error:.3e} V")
    return {"constant_current_charge_check": "passed", "constant_current_charge_max_error_V": max_error}


def main():
    output_dir = Path(__file__).resolve().parent / "results"
    output_dir.mkdir(exist_ok=True)
    metrics = current_source_example(output_dir)
    metrics.update(rc_example(output_dir))
    metrics.update(constant_current_charge_example(output_dir))
    metrics["reference_method"] = "Analytical formulas and backward Euler recurrence; MATLAB not run."
    (output_dir / "validation_summary.json").write_text(
        json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    print(f"Validation evidence saved in {output_dir}")


if __name__ == "__main__":
    main()
