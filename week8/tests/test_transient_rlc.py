"""Zhihang Chen -- three-stage RLC transient analysis.

Run from any directory: python tests/test_transient_rlc.py
Writes full-resolution BE/TR waveform CSVs for a constant 5 V source.
"""

import csv
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit_simulator"))
from capacitor import Capacitor
from inductor import Inductor
from linear_solvers import backward_euler, trapezoidal_rule
from mna_assembler import MNAAssembler
from resistor import Resistor
from voltage_source import VoltageSource


# Positive inductor current flows from node_plus to node_minus.
COMPONENTS = [
    ("R1", "in", "a1", 10.0),
    ("L1", "a1", "n1", 10e-3),
    ("C1", "n1", "gnd", 1e-6),
    ("R2", "n1", "a2", 15.0),
    ("L2", "a2", "n2", 22e-3),
    ("C2", "n2", "gnd", 2.2e-6),
    ("R3", "n2", "a3", 10.0),
    ("L3", "a3", "n3", 15e-3),
    ("C3", "n3", "gnd", 1.5e-6),
    ("R4", "n2", "gnd", 330.0),
    ("R5", "n3", "gnd", 1000.0),
]
TIME_STEPS_US = (1, 5, 10)
STOP_TIME_S = 30e-3
INPUT_V = 5.0
OUTPUT_DIR = Path(__file__).resolve().parent / "rlc_results"


def build_circuit():
    circuit = MNAAssembler()
    circuit.add_element(VoltageSource("Vin", "in", "gnd", INPUT_V))
    classes = {"R": Resistor, "L": Inductor, "C": Capacitor}
    for name, plus, minus, value in COMPONENTS:
        circuit.add_element(classes[name[0]](name, plus, minus, value))
    return circuit


def main():
    # 1. Build and assemble the circuit.
    circuit = build_circuit()
    G, M, b = circuit.assemble_system()
    node_indices = [circuit.node_map.get_index(f"n{i}") for i in range(1, 4)]
    inductor_indices = [circuit.node_map.get_branch_index(f"L{i}") for i in range(1, 4)]
    observed_indices = node_indices + inductor_indices
    labels = [f"V_n{i}_V" for i in range(1, 4)] + [f"I_L{i}_A" for i in range(1, 4)]

    # 2. Start with zero capacitor voltages and zero inductor currents.
    # Vin is already 5 V at t=0. L1 initially carries no current, so R1 has
    # no voltage drop: V(in)=V(a1)=5 V. Other node voltages and Vin current
    # are zero. The entire MNA vector must NOT be set to zero with a DC source.
    x0 = np.zeros(b.size)
    x0[circuit.node_map.get_index("in")] = INPUT_V
    x0[circuit.node_map.get_index("a1")] = INPUT_V

    def source_at(time):
        return b

    # 3. Run both existing solvers at every selected fixed step.
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for step_us in TIME_STEPS_US:
        h = step_us * 1e-6
        solutions = {}
        times = None
        for method, solve in (("BE", backward_euler), ("TR", trapezoidal_rule)):
            time, states = solve(G, M, source_at, x0.copy(), h, STOP_TIME_S)
            if times is not None and not np.array_equal(times, time):
                raise ValueError("BE and TR must use the same time grid.")
            if not np.array_equal(states[0], x0):
                raise ValueError("Solver did not retain the supplied initial state.")
            times = time
            solutions[method] = states[:, observed_indices]

        # 4. One CSV per step: same-time BE/TR pairs, full resolution, SI units.
        headers = ["time_s", "Vin_V"]
        for label in labels:
            headers.extend((f"BE_{label}", f"TR_{label}"))
        columns = [times, np.full(times.shape, INPUT_V)]
        for column in range(len(labels)):
            columns.extend((solutions["BE"][:, column], solutions["TR"][:, column]))
        data = np.column_stack(columns)
        path = OUTPUT_DIR / f"rlc_waveforms_h_{step_us}us.csv"
        with path.open("w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerow(headers)
            writer.writerows(data)
        print(f"h={step_us:2d} us: BE + TR, {len(times):5d} samples each -> {path.name}")

    print(f"Saved waveform CSVs to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
