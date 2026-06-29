#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Assemble the final ternary complex PDB:
    6UJN (protein) + Meropenem (first ligand) + Vancomycin (second ligand)
"""

from pathlib import Path

WORK = Path(r"D:\6UJN_Meropenem_Vancomycin")
OUTPUT = WORK / "output"

receptor = OUTPUT / "6UJN_Meropenem_rigid.pdbqt"
vancomycin = OUTPUT / "Vancomycin_out.pdbqt"
out_pdb = OUTPUT / "6UJN_Meropenem_Vancomycin.pdb"


def atom_lines(path, model=1):
    """Extract ATOM/HETATM lines from the first (or specified) model."""
    lines = []
    current_model = 0
    in_model = False
    has_models = False
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("MODEL"):
                has_models = True
                current_model += 1
                in_model = (current_model == model)
                continue
            if line.startswith("ENDMDL"):
                in_model = False
                continue
            if has_models and not in_model:
                continue
            if line.startswith("ATOM") or line.startswith("HETATM"):
                lines.append(line.rstrip("\n"))
    return lines


def renumber_atoms(lines):
    """Renumber ATOM/HETATM serial numbers sequentially."""
    out = []
    for i, line in enumerate(lines, 1):
        # PDB serial is columns 7-11 (5 digits)
        serial = min(i, 99999)
        new_line = line[:6] + f"{serial:5d}" + line[11:]
        out.append(new_line)
    return out


def main():
    rec_lines = atom_lines(receptor)
    # Vancomycin should be the second ligand; rename chain/residue to avoid confusion
    van_lines_raw = atom_lines(vancomycin)
    van_lines = []
    for line in van_lines_raw:
        # Set residue name to VAN and chain B, resseq 1
        new_line = line[:17] + "VAN" + " B" + "   1" + line[26:]
        van_lines.append(new_line)

    all_lines = rec_lines + van_lines
    all_lines = renumber_atoms(all_lines)

    with open(out_pdb, "w", encoding="utf-8") as f:
        f.write("REMARK   6UJN + Meropenem + Vancomycin ternary complex\n")
        f.write("REMARK   Meropenem-first sequential docking result\n")
        f.write("REMARK   Vancomycin best affinity: -10.144 kcal/mol\n")
        f.write("TER\n".join([""] * 0))
        for line in all_lines:
            f.write(line + "\n")
        f.write("TER\n")
        f.write("END\n")

    print(f"Ternary complex PDB saved: {out_pdb}")
    print(f"Total atoms: {len(all_lines)}")


if __name__ == "__main__":
    main()
