'''Zhihang Chen'''
'''Compare Backward Euler and the trapezoidal rule against the RC analytical solution.'''

import csv
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit_simulator"))
from capacitor import Capacitor
from linear_solvers import backward_euler, trapezoidal_rule
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


def main():
    # 1. Build the RC circuit.
    R = 1000.0
    C = 1e-6
    Vin = 5.0
    v0 = 0.0

    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vin", "in", "gnd", Vin))
    circuit.add_element(Resistor("R1", "in", "out", R))
    circuit.add_element(Capacitor("C1", "out", "gnd", C))
    G, M, b = circuit.assemble_system()

    node_in = circuit.node_map.get_index("in")
    node_out = circuit.node_map.get_index("out")
    source_current = circuit.node_map.get_branch_index("Vin")
    x0 = np.zeros(b.size)
    x0[node_in] = Vin
    x0[node_out] = v0
    x0[source_current] = -(Vin - v0) / R

    # 2. Choose the simulation settings shared by both methods.
    tau = R * C
    tstop = 5 * tau
    comparison_time = 2e-3
    h_values = [tau / 10, tau / 20, tau / 40, tau / 80, tau / 160]
    methods = [("BE", backward_euler), ("TR", trapezoidal_rule)]

    # 3. Compare each numerical solution with the analytical voltage.
    output_file = Path(__file__).with_name("test_transient_summary.csv")
    waveform_file = Path(__file__).with_name("test_transient_waveforms.csv")
    summary_rows = []
    waveform_rows = []

    for h in h_values:
        sample_index = round(comparison_time / h)
        for method, solve in methods:
            t, X = solve(G, M, lambda time: b, x0, h, tstop)
            if sample_index >= len(t) or not np.isclose(
                t[sample_index], comparison_time, rtol=1e-10, atol=0
            ):
                raise ValueError(f"h={h:g} s does not sample t={comparison_time:g} s.")

            v_numerical = X[:, node_out]
            v_analytical = Vin + (v0 - Vin) * np.exp(-t / tau)

            abs_error = np.abs(v_numerical - v_analytical)
            max_error = np.max(abs_error)
            rms_error = np.sqrt(np.mean(abs_error ** 2))
            error_2ms = abs_error[sample_index]
            summary_rows.append([h, method, max_error, rms_error, error_2ms])

            for i in range(len(t)):
                waveform_rows.append([h, method, t[i], v_analytical[i], v_numerical[i], abs_error[i]])

    with output_file.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["h_s", "method", "max_error_V", "rms_error_V", "2ms_error_V"])
        writer.writerows(summary_rows)

    with waveform_file.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow([
            "h_s", "method", "time_s", "analytical_voltage_V", "numerical_voltage_V", "abs_error_V"
        ])
        writer.writerows(waveform_rows)

    print(f"Summary saved to {output_file}")
    print(f"Waveform data saved to {waveform_file}")


if __name__ == "__main__":
    main()
