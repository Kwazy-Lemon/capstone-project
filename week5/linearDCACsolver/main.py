"""Run circuit.txt by default, or: python3 main.py ac_rc.txt"""
import argparse
import cmath
import math
from pathlib import Path
from netlist_reader import netlist_reader
from netlist_parser import parse_analyses
from MNA import mna
from ac_analysis import solve_ac


def run_dc(elements):
    G, b, x, unknowns = mna(elements)

    print("G matrix:")
    for row in G:
        print(row)

    print("b vector:", b)

    print("Solution:")
    for name, value in zip(unknowns, x):
        print(f"{name} = {value:.6g}")


def run_ac(elements, f_hz, ac_sources):
    A, b_ac, x_ac, unknowns = solve_ac(elements, f_hz, ac_sources)
    print(f"AC at {f_hz:g} Hz:")
    print("A(jw) matrix:")
    for row in A:
        print(row)
    print("b(jw) vector:", b_ac)
    print("AC solution:")
    for name, value in zip(unknowns, x_ac):
        unit = "V" if name.startswith("V(") else "A"
        phase = f"{math.degrees(cmath.phase(value)):.3f} deg" if value else "undefined"
        print(f"{name} = {value:.6g} {unit}; magnitude = {abs(value):.6g} {unit}; phase = {phase}")


def main():
    parser = argparse.ArgumentParser(description="Linear DC and single-frequency AC analysis")
    parser.add_argument("netlist", nargs="?", default=Path(__file__).with_name("circuit.txt"))
    entries = netlist_reader(parser.parse_args().netlist)
    elements, dc_requested, f_hz, ac_sources = parse_analyses(entries)
    print("DC requested:", dc_requested)
    for element in elements:
        print(element)

    failed = False
    if dc_requested:
        try:
            run_dc(elements)
        except (ValueError, ZeroDivisionError) as err:
            print(f"DC error: {err}")
            failed = True
    if f_hz is not None:
        try:
            run_ac(elements, f_hz, ac_sources)
        except ValueError as err:
            print(f"AC error: {err}")
            failed = True
    if not dc_requested and f_hz is None:
        print("No analysis requested. Add .dc or .ac F_HZ before .end.")
    return int(failed)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as err:
        raise SystemExit(f"Input error: {err}")

