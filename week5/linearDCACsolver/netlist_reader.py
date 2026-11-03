from pathlib import Path

def netlist_reader(file_path):
    lines = Path(file_path).read_text(
        encoding="utf-8"
    ).splitlines()

    entries = []                 
    for line in lines[1:]:
        line = line.strip()

        if not line:
            continue

        if line.lower() == ".end":
            break

        entries.append(line)     

    return entries              
