"""Convert this project's netlist lines into elements and analysis requests."""
import cmath
import math

from elements import (Resistor, Capacitor, Inductor, VoltageSource, CurrentSource,
                      VCCS, VCVS, CCCS, CCVS)


def number(token):
    """Use the same decimal/scientific notation as the original DC parser."""
    value = float(token)
    if not math.isfinite(value):
        raise ValueError(f"Expected a finite number, got {token!r}")
    return value


def source_values(fields):
    """Read DC value and/or AC magnitude [phase_deg] on an independent source."""
    dc_value, ac_phasor = 0.0, None
    seen = set()
    i = 0
    while i < len(fields):
        mode = fields[i].upper()
        if mode not in ("DC", "AC") or mode in seen or i + 1 >= len(fields):
            raise ValueError("Expected DC value and/or AC magnitude [phase_deg]")
        seen.add(mode)
        value = number(fields[i + 1])
        i += 2
        if mode == "DC":
            dc_value = value
        else:
            if value < 0:
                raise ValueError("AC magnitude must be nonnegative")
            phase = 0.0
            if i < len(fields) and fields[i].upper() not in ("DC", "AC"):
                phase = number(fields[i])
                i += 1
            ac_phasor = cmath.rect(value, math.radians(phase))
    return dc_value, ac_phasor


def parse_analyses(entries):
    """Return elements, DC flag, AC frequency in Hz (or None), AC source phasors."""
    elements = []
    dc_requested = False
    f_hz = None
    ac_sources = {}
    element_names = set()

    for entry in entries:
        parts = entry.split()
        if not parts:
            continue
        name = parts[0].upper()
        try:
            if name == ".DC":
                if len(parts) != 1:
                    raise ValueError("Use .dc by itself for a DC operating point")
                dc_requested = True
                continue
            if name == ".AC":
                if len(parts) != 2 or f_hz is not None:
                    raise ValueError("Use one .ac F_HZ line")
                f_hz = number(parts[1])
                if f_hz <= 0:
                    raise ValueError("AC frequency must be positive")
                continue

            if name in element_names:
                raise ValueError(f"Duplicate element name: {parts[0]}")
            # Test the long names before V for voltage source or C for capacitor.
            if name.startswith(("VCCS", "VCVS")) and len(parts) == 6:
                cls = VCCS if name.startswith("VCCS") else VCVS
                element = cls(parts[0], parts[1], parts[2], parts[3], parts[4], number(parts[5]))
            elif name.startswith(("CCCS", "CCVS")) and len(parts) == 5:
                cls = CCCS if name.startswith("CCCS") else CCVS
                element = cls(parts[0], parts[1], parts[2], parts[3], number(parts[4]))
            elif name.startswith(("VCCS", "VCVS", "CCCS", "CCVS")):
                raise ValueError("Invalid controlled-source format")
            elif name.startswith(("R", "C", "L")) and len(parts) == 4:
                cls = {"R": Resistor, "C": Capacitor, "L": Inductor}[name[0]]
                element = cls(parts[0], parts[1], parts[2], number(parts[3]))
            elif name.startswith(("V", "I")) and len(parts) >= 5:
                dc_value, ac_phasor = source_values(parts[3:])
                cls = VoltageSource if name.startswith("V") else CurrentSource
                element = cls(parts[0], parts[1], parts[2], dc_value)
                if ac_phasor is not None:
                    ac_sources[name] = ac_phasor
            else:
                raise ValueError("Unsupported element or invalid format")
        except ValueError as error:
            raise ValueError(f"{entry}: {error}") from error
        element_names.add(name)
        elements.append(element)

    return elements, dc_requested, f_hz, ac_sources


def netlist_parser(entries):
    """Original two-result interface for existing DC callers."""
    elements, dc_requested, _, _ = parse_analyses(entries)
    return elements, dc_requested
