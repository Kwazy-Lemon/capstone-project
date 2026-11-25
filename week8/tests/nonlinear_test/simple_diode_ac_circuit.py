"""Small-signal AC test of a DC-biased diode with a parallel capacitor."""

import csv
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "circuit_simulator"))

from capacitor import Capacitor
from mna_assembler import MNAAssembler
from nonlinear.diode import Diode
from nonlinear.nonlinear_solver import NonlinearSolver
from resistor import Resistor
from voltage_source import VoltageSource


def main(plot=False):
    supply_voltage = 5.0
    ac_amplitude = 1e-3
    resistance = 1000.0
    capacitance = 1e-6
    saturation_current = 1e-9
    ideality_factor = 1.0
    thermal_voltage = 0.0258649529170901  # Match LTspice at 27 degrees C.
    frequencies_hz = np.logspace(1, 6, 101)

    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("V1", "in", "0", supply_voltage))
    circuit.add_element(Resistor("R1", "in", "out", resistance))
    circuit.add_element(
        Diode("D1", "out", "0", Is=saturation_current,
              n=ideality_factor, Vt=thermal_voltage)
    )
    circuit.add_element(Capacitor("C1", "out", "0", capacitance))
    G, _, _ = circuit.assemble_system()
    node_in = circuit.node_map.get_index("in")
    node_out = circuit.node_map.get_index("out")
    source_index = circuit.node_map.get_branch_index("V1")

    # DC component values set the bias. This separate vector sets the AC input.
    b_ac = np.zeros(G.shape[0], dtype=complex)
    b_ac[source_index] = ac_amplitude
    x0 = np.zeros(G.shape[0])
    x0[node_in] = supply_voltage
    x0[node_out] = 0.4
    x0[source_index] = -(supply_voltage - 0.4) / resistance

    solver = NonlinearSolver(circuit)
    x_dc, x_ac = solver.ac_sweep(
        b_ac, frequencies_hz, x0=x0, tol=1e-9, residual_tol=1e-12
    )

    # Independent scalar DC reference, without MNA or the Newton solver.
    def current_difference(voltage):
        return saturation_current * np.expm1(
            voltage / (ideality_factor * thermal_voltage)
        ) - (supply_voltage - voltage) / resistance

    expected_dc_voltage = brentq(current_difference, 0, supply_voltage, xtol=1e-14)
    expected_dc_current = -(supply_voltage - expected_dc_voltage) / resistance
    np.testing.assert_allclose(x_dc[[node_in, node_out]],
                               [supply_voltage, expected_dc_voltage],
                               rtol=0, atol=1e-9)
    np.testing.assert_allclose(x_dc[source_index], expected_dc_current,
                               rtol=0, atol=1e-12)

    # KCL gives Vout/Vin = (1/R) / (1/R + gd + j*omega*C).
    gd = saturation_current / (ideality_factor * thermal_voltage) * np.exp(
        expected_dc_voltage / (ideality_factor * thermal_voltage)
    )
    omega = 2 * np.pi * frequencies_hz
    expected_output = ac_amplitude / (
        1 + resistance * gd + 1j * omega * resistance * capacitance
    )
    expected_source = -(ac_amplitude - expected_output) / resistance
    np.testing.assert_allclose(x_ac[:, node_in], ac_amplitude, rtol=1e-10, atol=1e-15)
    np.testing.assert_allclose(x_ac[:, node_out], expected_output, rtol=1e-9, atol=1e-15)
    np.testing.assert_allclose(x_ac[:, source_index], expected_source, rtol=1e-9, atol=1e-15)

    simulated = {
        "V(in)": x_ac[:, node_in],
        "V(out)": x_ac[:, node_out],
        "I(V1)": x_ac[:, source_index],
        "I(R1)": (x_ac[:, node_in] - x_ac[:, node_out]) / resistance,
        "I(D1)": gd * x_ac[:, node_out],
        "I(C1)": 1j * omega * capacitance * x_ac[:, node_out],
    }
    reference_path = Path(__file__).with_name("simple_diode_ac_ltspice_reference.csv")
    with reference_path.open(newline="", encoding="utf-8") as reference_file:
        rows = list(csv.DictReader(reference_file))
    np.testing.assert_allclose(frequencies_hz, [float(row["frequency_hz"]) for row in rows],
                               rtol=1e-10, atol=0)
    for name, values in simulated.items():
        expected = np.array([
            complex(float(row[name + "_real"]), float(row[name + "_imag"]))
            for row in rows
        ])
        np.testing.assert_allclose(values, expected, rtol=1e-6, atol=1e-14,
                                   err_msg=f"LTspice AC comparison for {name}")

    gain = x_ac[:, node_out] / x_ac[:, node_in]
    pole_hz = (1 / resistance + gd) / (2 * np.pi * capacitance)
    print("PASS: simple diode small-signal AC analysis")
    print(f"DC diode voltage: {x_dc[node_out]:.9f} V")
    print(f"Diode differential conductance: {gd:.9f} S")
    print(f"Analytical pole frequency: {pole_hz:.3f} Hz")
    print(f"{'Frequency (Hz)':>15} {'|Vout| (uV)':>15} {'Gain (dB)':>12} {'Phase (deg)':>14}")
    for index in range(0, len(frequencies_hz), 20):
        print(f"{frequencies_hz[index]:15.1f} {abs(x_ac[index, node_out]) * 1e6:15.6f}"
              f" {20 * np.log10(abs(gain[index])):12.6f} {np.angle(gain[index], deg=True):14.6f}")
    print("All 101 frequencies passed the analytical and saved LTspice comparisons.")
    if plot:
        import matplotlib.pyplot as plt

        gain_db = 20 * np.log10(np.abs(gain))
        phase_deg = np.rad2deg(np.unwrap(np.angle(gain)))
        fig, (magnitude_ax, phase_ax) = plt.subplots(
            2, 1, figsize=(8, 6.5), sharex=True
        )
        magnitude_ax.semilogx(frequencies_hz, gain_db, linewidth=2)
        magnitude_ax.set_ylabel("Gain (dB)")
        magnitude_ax.set_title("Simple diode Bode plot")
        phase_ax.semilogx(frequencies_hz, phase_deg, linewidth=2)
        phase_ax.set_xlabel("Frequency (Hz)")
        phase_ax.set_ylabel("Phase (degrees)")
        for ax in (magnitude_ax, phase_ax):
            ax.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        plot_path = Path(__file__).with_name("simple_diode_ac_response.png")
        fig.savefig(plot_path, dpi=180)
        print(f"Plot saved to {plot_path}")
        plt.show()

    return x_dc, x_ac


if __name__ == "__main__":
    main(plot=True)
