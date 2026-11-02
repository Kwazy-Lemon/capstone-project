'''Zhihang Chen'''
'''Compare the RC analytical solution and Backward Euler for several h values.'''

import csv
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit_simulator"))
from capacitor import Capacitor
from linear_solvers import backward_euler
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

    # 2. Choose the end time and several fixed time steps.
    tau = R * C
    tstop = 5 * tau
    h_values = [tau / 10, tau / 20, tau / 40, tau / 80, tau / 160]

    # 3. Compare both solutions for each h and record their sampled values.
    output_file = Path(__file__).with_name("test_transient_summary.csv")
    waveform_file = Path(__file__).with_name("test_transient_waveforms.csv")
    waveform_rows = []
    with output_file.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["h_s", "max_error_V", "rms_error_V"])

        for h in h_values:
            t, X = backward_euler(G, M, lambda time: b, x0, h, tstop)
            v_be = X[:, node_out]
            v_analytical = Vin + (v0 - Vin) * np.exp(-t / tau)

            error = v_be - v_analytical
            max_error = np.max(np.abs(error))
            rms_error = np.sqrt(np.mean(error ** 2))
            writer.writerow([h, max_error, rms_error])

            for i in range(len(t)):
                waveform_rows.append([h, t[i], v_analytical[i], v_be[i], abs(error[i])])

    with waveform_file.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["h_s", "time_s", "analytical_voltage_V", "backward_euler_voltage_V", "abs_error_V"])
        writer.writerows(waveform_rows)

    print(f"Summary saved to {output_file}")
    print(f"Waveform data saved to {waveform_file}")


if __name__ == "__main__":
    main()
