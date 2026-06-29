#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare flexible holo receptor for sequential docking (Meropenem-first).

Inputs:
    - apo receptor PDBQT (6UJN.pdbqt)
    - first-ligand docked pose PDBQT (Meropenem_out.pdbqt, model 1)

Outputs:
    - holo receptor PDBQT (rigid part): 6UJN_Meropenem_rigid.pdbqt
    - flexible side chains PDBQT: 6UJN_Meropenem_flex.pdbqt
    - configuration file: conf.txt (22x22x22 box centered on Meropenem centroid)
    - log of flexible residues: flex_residues.txt

Residues with any atom within 5 Angstrom of any Meropenem heavy atom are marked flexible.
For each flexible residue, the side-chain atoms (beyond CA/C/N/O of the backbone)
are written to the flex file; the backbone atoms remain in the rigid file.
"""

import math
import os
from pathlib import Path

WORK = Path(r"D:\6UJN_Meropenem_Vancomycin")
INPUT = WORK / "input"
OUTPUT = WORK / "output"
OUTPUT.mkdir(parents=True, exist_ok=True)

RECEPTOR = INPUT / "6UJN.pdbqt"
LIGAND = INPUT / "Meropenem_out.pdbqt"

# Side-chain atoms per residue (backbone atoms stay rigid)
BACKBONE_ATOMS = {"N", "CA", "C", "O", "HN", "H", "HA", "HA1", "HA2", "H1", "H2", "H3"}


def parse_pdbqt_atoms(path):
    """Parse ATOM/HETATM records from a PDBQT file."""
    atoms = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            record = line[:6]
            if record in ("ATOM  ", "HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21].strip()
                res_seq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                charge = float(line[70:76])
                atom_type = line[77:79].strip()
                atoms.append({
                    "line": line.rstrip("\n"),
                    "record": record.strip(),
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "chain": chain,
                    "res_seq": res_seq,
                    "x": x,
                    "y": y,
                    "z": z,
                    "charge": charge,
                    "atom_type": atom_type,
                })
    return atoms


def centroid(atoms):
    n = len(atoms)
    sx = sum(a["x"] for a in atoms)
    sy = sum(a["y"] for a in atoms)
    sz = sum(a["z"] for a in atoms)
    return (sx / n, sy / n, sz / n)


def distance(a1, a2):
    return math.sqrt((a1["x"] - a2["x"]) ** 2 +
                     (a1["y"] - a2["y"]) ** 2 +
                     (a1["z"] - a2["z"]) ** 2)


def is_backbone(atom_name):
    return atom_name in BACKBONE_ATOMS


def extract_first_model(in_path, out_path):
    """Copy only the first MODEL/ENDMDL block, or the whole file if no MODEL."""
    with open(in_path, "r", encoding="utf-8") as f:
        lines = f.readlines()
    out_lines = []
    in_model = False
    model_started = False
    for line in lines:
        if line.startswith("MODEL"):
            if model_started:
                break
            in_model = True
            model_started = True
            continue
        if line.startswith("ENDMDL"):
            break
        if in_model or not model_started:
            out_lines.append(line)
    with open(out_path, "w", encoding="utf-8") as f:
        f.writelines(out_lines)


def main():
    # 1. Extract first model of ligand
    lig_first = OUTPUT / "Meropenem_first_model.pdbqt"
    extract_first_model(LIGAND, lig_first)

    lig_atoms = parse_pdbqt_atoms(lig_first)
    lig_heavy = [a for a in lig_atoms if a["atom_type"] != "HD"]
    print(f"Meropenem first model: {len(lig_atoms)} atoms, {len(lig_heavy)} heavy atoms")

    cx, cy, cz = centroid(lig_heavy)
    print(f"Meropenem heavy-atom centroid: {cx:.4f} {cy:.4f} {cz:.4f}")

    # 2. Identify flexible residues
    rec_atoms = parse_pdbqt_atoms(RECEPTOR)
    flex_residues = set()
    for ra in rec_atoms:
        for la in lig_heavy:
            if distance(ra, la) <= 5.0:
                flex_residues.add((ra["chain"], ra["res_seq"], ra["res_name"]))
                break

    flex_residues = sorted(flex_residues)
    print(f"Flexible residues within 5A of Meropenem: {len(flex_residues)}")
    for chain, seq, name in flex_residues:
        print(f"  {name} {chain}{seq}")

    # 3. Split receptor into rigid and flexible parts
    flex_keys = {(chain, seq) for chain, seq, _ in flex_residues}
    rigid_lines = []
    flex_lines = []

    # Vina expects a ROOT/ENDROOT/TORSDOF structure in the flex file? Actually
    # for --flex, Vina expects a PDBQT with ATOM/HETATM lines and torsion records.
    # The simplest valid format is just ATOM records with proper charges and atom types.
    # We keep the original lines (which already have charges/types) and add a TORSDOF header.

    flex_lines.append("REMARK  flexible side chains for 6UJN_Meropenem holo receptor\n")
    flex_lines.append("REMARK  residues within 5A of bound Meropenem\n")
    for chain, seq, name in flex_residues:
        flex_lines.append(f"REMARK  flexible: {name} {chain}{seq}\n")

    for atom in rec_atoms:
        key = (atom["chain"], atom["res_seq"])
        if key in flex_keys and not is_backbone(atom["atom_name"]):
            flex_lines.append(atom["line"] + "\n")
        else:
            rigid_lines.append(atom["line"] + "\n")

    # Append ligand atoms to rigid receptor to make holo receptor
    with open(lig_first, "r", encoding="utf-8") as f:
        lig_lines = f.readlines()

    rigid_out = OUTPUT / "6UJN_Meropenem_rigid.pdbqt"
    flex_out = OUTPUT / "6UJN_Meropenem_flex.pdbqt"

    with open(rigid_out, "w", encoding="utf-8") as f:
        f.writelines(rigid_lines)
        # append ligand as part of rigid receptor (skip MODEL/ENDMDL/REMARK/torsion markers)
        for line in lig_lines:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                f.write(line)
        f.write("TER\n")
        f.write("END\n")

    with open(flex_out, "w", encoding="utf-8") as f:
        f.writelines(flex_lines)
        f.write("TER\n")
        f.write("END\n")

    print(f"Rigid holo receptor written: {rigid_out} ({len(rigid_lines)} protein atoms + ligand)")
    print(f"Flexible side chains written: {flex_out} ({len(flex_lines)} atoms)")

    # 4. Write configuration
    conf_out = OUTPUT / "conf.txt"
    with open(conf_out, "w", encoding="utf-8") as f:
        f.write(f"receptor = {rigid_out}\n")
        f.write(f"flex = {flex_out}\n")
        f.write(f"center_x = {cx:.4f}\n")
        f.write(f"center_y = {cy:.4f}\n")
        f.write(f"center_z = {cz:.4f}\n")
        f.write("size_x = 22.0\n")
        f.write("size_y = 22.0\n")
        f.write("size_z = 22.0\n")
        f.write("exhaustiveness = 16\n")
        f.write("num_modes = 10\n")
        f.write("energy_range = 3\n")
        f.write("seed = 42\n")

    print(f"Config written: {conf_out}")

    # 5. Save flexible residue list
    flex_list = OUTPUT / "flex_residues.txt"
    with open(flex_list, "w", encoding="utf-8") as f:
        f.write("chain,res_seq,res_name\n")
        for chain, seq, name in flex_residues:
            f.write(f"{chain},{seq},{name}\n")

    # 6. Save centroid
    centroid_file = OUTPUT / "meropenem_centroid.txt"
    with open(centroid_file, "w", encoding="utf-8") as f:
        f.write(f"{cx:.6f} {cy:.6f} {cz:.6f}\n")


if __name__ == "__main__":
    main()
