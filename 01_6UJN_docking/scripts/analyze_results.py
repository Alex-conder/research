#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Analyze sequential docking results for 6UJN Meropenem-first vs Vancomycin-first.

Outputs:
    - CSV of binding modes from the new Meropenem-first Vancomycin docking
    - RMSD of Vancomycin in ternary complex vs apo reference
    - Hydrogen bond counts for key residues (Tyr, Lys, etc.)
    - Cross-comparison report (Markdown table)
"""

import csv
import math
import os
import re
from pathlib import Path

WORK = Path(r"D:\6UJN_Meropenem_Vancomycin")
INPUT = WORK / "input"
OUTPUT = WORK / "output"

# Reference apo Vancomycin pose (single docking into 6UJN)
APO_VAN_DOCKED = Path(r"D:\Tools\Docking\docking_project\data\receptors\6UJN\Vancomycin_out.pdbqt")

# New result
NEW_VAN_DOCKED = OUTPUT / "Vancomycin_out.pdbqt"
NEW_LOG = OUTPUT / "log_rigid.txt"

# Existing ternary complex (Vancomycin-first)
VAN_FIRST_MERO_DOCKED = Path(r"D:\Tools\Docking\docking_project\data\receptors\6UJN_Vancomycin\Meropenem_out.pdbqt")

# Receptors
APO_RECEPTOR = INPUT / "6UJN.pdbqt"
MERO_HOLO_RECEPTOR = OUTPUT / "6UJN_Meropenem_rigid.pdbqt"
VAN_HOLO_RECEPTOR = Path(r"D:\Tools\Docking\docking_project\data\receptors\6UJN_Vancomycin\6UJN_Vancomycin.pdbqt")

# Key residues to monitor (D-Ala-D-Ala recognition: Tyr, Lys, etc.)
KEY_RES_NAMES = {"TYR", "LYS", "ASN", "GLN", "SER", "THR", "ARG", "HIS", "ASP", "GLU", "TRP"}

# H-bond geometric criteria
DONOR_ACCEP_MAX = 3.5      # D...A distance (Angstrom)
HYDROGEN_ACCEP_MAX = 2.5   # H...A distance (Angstrom)
ANGLE_MIN = 120.0          # D-H...A angle (degrees)
BOND_CUTOFF = 1.2          # X-H bond length for associating H to heavy atom


def parse_pdbqt_atoms(path, model=1):
    """Parse ATOM/HETATM from a PDBQT file, selecting the requested MODEL."""
    atoms = []
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
            record = line[:6]
            if record in ("ATOM  ", "HETATM"):
                atoms.append({
                    "line": line.rstrip("\n"),
                    "atom_name": line[12:16].strip(),
                    "res_name": line[17:20].strip(),
                    "chain": line[21].strip(),
                    "res_seq": int(line[22:26]),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "atom_type": line[77:79].strip(),
                })
    return atoms


def dist(a, b):
    return math.sqrt((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2)


def angle(p_vertex, p1, p2):
    """Angle p1 - p_vertex - p2 in degrees."""
    v1 = (p1["x"] - p_vertex["x"], p1["y"] - p_vertex["y"], p1["z"] - p_vertex["z"])
    v2 = (p2["x"] - p_vertex["x"], p2["y"] - p_vertex["y"], p2["z"] - p_vertex["z"])
    d1 = math.sqrt(sum(c ** 2 for c in v1))
    d2 = math.sqrt(sum(c ** 2 for c in v2))
    if d1 == 0 or d2 == 0:
        return 0.0
    dot = sum(v1[i] * v2[i] for i in range(3))
    cos = max(-1.0, min(1.0, dot / (d1 * d2)))
    return math.degrees(math.acos(cos))


def is_heavy(atom):
    t = atom["atom_type"]
    return t not in ("H", "HD", "HS", "D")


def is_potential_donor_heavy(atom):
    """N or O that can bear an H (donor)."""
    t = atom["atom_type"]
    return t.startswith("N") or t.startswith("O")


def is_potential_acceptor_heavy(atom):
    """N or O that can accept an H-bond."""
    t = atom["atom_type"]
    return t.startswith("N") or t.startswith("O") or t == "S" or t == "SA" or t == "Cl"


def assign_hydrogens(heavy_atoms, all_atoms):
    """Map each heavy atom to its covalently bonded hydrogens (distance < 1.2 A)."""
    hydrogens = [a for a in all_atoms if not is_heavy(a)]
    mapping = {id(h): None for h in hydrogens}
    heavy_h = {id(h): [] for h in heavy_atoms}
    for ha in heavy_atoms:
        for h in hydrogens:
            if dist(ha, h) < BOND_CUTOFF:
                heavy_h[id(ha)].append(h)
                mapping[id(h)] = ha
    return heavy_h, mapping


def count_hbonds(rec_atoms, lig_atoms):
    """Count hydrogen bonds between receptor and ligand heavy atoms."""
    rec_heavy = [a for a in rec_atoms if is_heavy(a)]
    lig_heavy = [a for a in lig_atoms if is_heavy(a)]
    rec_h_map, _ = assign_hydrogens(rec_heavy, rec_atoms)
    lig_h_map, _ = assign_hydrogens(lig_heavy, lig_atoms)

    count = 0
    details = []

    # Receptor donor -> Ligand acceptor
    for ra in rec_heavy:
        if not is_potential_donor_heavy(ra):
            continue
        for la in lig_heavy:
            if not is_potential_acceptor_heavy(la):
                continue
            if dist(ra, la) > DONOR_ACCEP_MAX:
                continue
            for h in rec_h_map.get(id(ra), []):
                if dist(h, la) <= HYDROGEN_ACCEP_MAX:
                    ang = angle(h, ra, la)
                    if ang >= ANGLE_MIN:
                        count += 1
                        details.append(("R", ra, "L", la, dist(ra, la), dist(h, la), ang))
                        break

    # Ligand donor -> Receptor acceptor
    for la in lig_heavy:
        if not is_potential_donor_heavy(la):
            continue
        for ra in rec_heavy:
            if not is_potential_acceptor_heavy(ra):
                continue
            if dist(la, ra) > DONOR_ACCEP_MAX:
                continue
            for h in lig_h_map.get(id(la), []):
                if dist(h, ra) <= HYDROGEN_ACCEP_MAX:
                    ang = angle(h, la, ra)
                    if ang >= ANGLE_MIN:
                        count += 1
                        details.append(("L", la, "R", ra, dist(la, ra), dist(h, ra), ang))
                        break

    return count, details


def key_residue_hbonds(receptor_atoms, ligand_atoms):
    """Return H-bond counts per key residue."""
    key_rec = [a for a in receptor_atoms if a["res_name"] in KEY_RES_NAMES]
    _, details = count_hbonds(key_rec, ligand_atoms)
    counts = {}
    for donor_flag, donor, acceptor_flag, acceptor, d, dh, ang in details:
        if donor_flag == "R":
            res = (donor["chain"], donor["res_seq"], donor["res_name"])
        else:
            res = (acceptor["chain"], acceptor["res_seq"], acceptor["res_name"])
        counts[res] = counts.get(res, 0) + 1
    return counts


def rmsd(atoms1, atoms2):
    """Calculate RMSD between two lists of atoms; match by atom name if counts differ."""
    if len(atoms1) != len(atoms2):
        m1 = {a["atom_name"]: a for a in atoms1}
        m2 = {a["atom_name"]: a for a in atoms2}
        common = set(m1.keys()) & set(m2.keys())
        if not common:
            return None
        atoms1 = [m1[n] for n in common]
        atoms2 = [m2[n] for n in common]
    n = len(atoms1)
    ss = sum(dist(atoms1[i], atoms2[i]) ** 2 for i in range(n))
    return math.sqrt(ss / n)


def parse_vina_log(log_path):
    """Parse Vina log for binding modes table."""
    modes = []
    in_table = False
    with open(log_path, "r", encoding="utf-8") as f:
        for line in f:
            if "mode |" in line and "affinity" in line:
                in_table = True
                continue
            if in_table:
                parts = line.strip().split()
                if len(parts) >= 4 and parts[0].isdigit():
                    try:
                        modes.append({
                            "mode": int(parts[0]),
                            "affinity": float(parts[1]),
                            "rmsd_lb": float(parts[2]),
                            "rmsd_ub": float(parts[3]),
                        })
                    except ValueError:
                        continue
                elif line.strip() == "":
                    break
    return modes


def extract_ligand_from_holo(holo_pdbqt, output_pdbqt):
    """Extract HETATM records (the first ligand, Meropenem) from a holo receptor PDBQT."""
    hetatms = []
    with open(holo_pdbqt, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("HETATM"):
                hetatms.append(line)
    with open(output_pdbqt, "w", encoding="utf-8") as f:
        f.write("MODEL 1\n")
        for line in hetatms:
            f.write(line)
        f.write("ENDMDL\n")
        f.write("END\n")


def load_csv_best(csv_path):
    """Return the best (lowest) binding affinity from a Vina CSV summary."""
    best = None
    with open(csv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                aff = float(row["Binding Affinity"])
            except (KeyError, ValueError):
                continue
            if best is None or aff < best:
                best = aff
    return best


def main():
    print("=" * 60)
    print("6UJN Meropenem-first sequential docking analysis")
    print("=" * 60)

    # Existing reference values
    e_mero_apo = load_csv_best(INPUT / "6UJN_Meropenem.csv")
    e_van_apo = load_csv_best(INPUT / "6UJN_Vancomycin.csv")
    e_van_first_mero = load_csv_best(INPUT / "6UJN_Vancomycin_Meropenem.csv")

    print(f"Reference (apo) Meropenem affinity: {e_mero_apo:.3f} kcal/mol")
    print(f"Reference (apo) Vancomycin affinity: {e_van_apo:.3f} kcal/mol")
    print(f"Vancomycin-first ternary Meropenem affinity: {e_van_first_mero:.3f} kcal/mol")

    if not NEW_LOG.exists():
        print(f"Log file not found: {NEW_LOG}. Docking may still be running.")
        return

    modes = parse_vina_log(NEW_LOG)
    print(f"\nParsed {len(modes)} binding modes from {NEW_LOG}")
    if modes:
        print(f"Best affinity (Meropenem-first ternary Vancomycin): {modes[0]['affinity']:.3f} kcal/mol")

    e_mero_first_van = modes[0]["affinity"] if modes else None

    # Save CSV
    csv_path = OUTPUT / "6UJN_Meropenem_Vancomycin.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Ligand", "Binding Affinity", "rmsd/ub", "rmsd/lb"])
        for m in modes:
            writer.writerow(["6UJN_Meropenem_Vancomycin", m["affinity"], m["rmsd_ub"], m["rmsd_lb"]])
    print(f"Saved CSV: {csv_path}")

    # RMSD of Vancomycin pose in new ternary vs apo reference
    rmsd_van = None
    if NEW_VAN_DOCKED.exists() and APO_VAN_DOCKED.exists():
        new_van = parse_pdbqt_atoms(NEW_VAN_DOCKED, model=1)
        ref_van = parse_pdbqt_atoms(APO_VAN_DOCKED, model=1)
        new_van_heavy = [a for a in new_van if is_heavy(a)]
        ref_van_heavy = [a for a in ref_van if is_heavy(a)]
        rmsd_van = rmsd(new_van_heavy, ref_van_heavy)
        print(f"Vancomycin heavy-atom RMSD (Meropenem-first ternary vs apo): {rmsd_van:.3f} A")
        with open(OUTPUT / "rmsd_vancomycin.txt", "w", encoding="utf-8") as f:
            f.write(f"{rmsd_van:.6f}\n")

    # RMSD of Meropenem pose in ternary (Meropenem-first) vs apo reference
    rmsd_mero = None
    if (INPUT / "Meropenem_out.pdbqt").exists() and VAN_FIRST_MERO_DOCKED.exists():
        mero_first = parse_pdbqt_atoms(INPUT / "Meropenem_out.pdbqt", model=1)
        mero_in_van_first = parse_pdbqt_atoms(VAN_FIRST_MERO_DOCKED, model=1)
        mero_first_heavy = [a for a in mero_first if is_heavy(a)]
        mero_van_first_heavy = [a for a in mero_in_van_first if is_heavy(a)]
        rmsd_mero = rmsd(mero_first_heavy, mero_van_first_heavy)
        print(f"Meropenem heavy-atom RMSD (Meropenem-first vs Vancomycin-first): {rmsd_mero:.3f} A")
        with open(OUTPUT / "rmsd_meropenem.txt", "w", encoding="utf-8") as f:
            f.write(f"{rmsd_mero:.6f}\n")

    # Hydrogen bond analysis
    # Use the full holo receptor (protein + first ligand) for context
    rec_full = parse_pdbqt_atoms(MERO_HOLO_RECEPTOR)
    rec_protein = [a for a in rec_full if a["res_name"]]  # protein atoms have residue names
    new_van = parse_pdbqt_atoms(NEW_VAN_DOCKED, model=1)

    hb_total, details = count_hbonds(rec_protein, new_van)
    key_hb = key_residue_hbonds(rec_protein, new_van)
    print(f"\nTotal protein-Vancomycin H-bonds in Meropenem-first complex: {hb_total}")
    print("Top key residue H-bond counts (Vancomycin):")
    for res, cnt in sorted(key_hb.items(), key=lambda x: -x[1])[:10]:
        print(f"  {res[2]} {res[0]}{res[1]}: {cnt}")

    # Meropenem H-bonds in the ternary complex
    mero_holo = OUTPUT / "Meropenem_from_holo.pdbqt"
    extract_ligand_from_holo(MERO_HOLO_RECEPTOR, mero_holo)
    mero_atoms = parse_pdbqt_atoms(mero_holo)
    hb_mero_total, details_mero = count_hbonds(rec_protein, mero_atoms)
    key_hb_mero = key_residue_hbonds(rec_protein, mero_atoms)
    print(f"\nTotal protein-Meropenem H-bonds in ternary: {hb_mero_total}")
    print("Top key residue H-bond counts (Meropenem):")
    for res, cnt in sorted(key_hb_mero.items(), key=lambda x: -x[1])[:10]:
        print(f"  {res[2]} {res[0]}{res[1]}: {cnt}")

    # For comparison: Vancomycin-first ternary - Meropenem H-bonds
    van_holo = VAN_HOLO_RECEPTOR
    van_first_rec_protein = parse_pdbqt_atoms(van_holo)
    van_first_rec_protein = [a for a in van_first_rec_protein if a["res_name"]]
    mero_in_van_first = parse_pdbqt_atoms(VAN_FIRST_MERO_DOCKED, model=1)
    hb_mero_van_first_total, _ = count_hbonds(van_first_rec_protein, mero_in_van_first)
    key_hb_mero_van_first = key_residue_hbonds(van_first_rec_protein, mero_in_van_first)
    print(f"\nTotal protein-Meropenem H-bonds in Vancomycin-first complex: {hb_mero_van_first_total}")
    print("Top key residue H-bond counts (Meropenem, Van-first):")
    for res, cnt in sorted(key_hb_mero_van_first.items(), key=lambda x: -x[1])[:10]:
        print(f"  {res[2]} {res[0]}{res[1]}: {cnt}")

    # Save H-bond summary
    with open(OUTPUT / "hbond_summary.txt", "w", encoding="utf-8") as f:
        f.write("6UJN sequential docking H-bond summary\n")
        f.write("=" * 60 + "\n\n")
        f.write("Meropenem-first ternary complex\n")
        f.write(f"  Protein-Vancomycin H-bonds: {hb_total}\n")
        f.write(f"  Protein-Meropenem H-bonds: {hb_mero_total}\n")
        f.write("\nVancomycin-first ternary complex\n")
        f.write(f"  Protein-Meropenem H-bonds: {hb_mero_van_first_total}\n")
        f.write("\nVancomycin key residues (Meropenem-first):\n")
        for res, cnt in sorted(key_hb.items(), key=lambda x: -x[1]):
            f.write(f"  {res[2]} {res[0]}{res[1]}: {cnt}\n")
        f.write("\nMeropenem key residues (Meropenem-first):\n")
        for res, cnt in sorted(key_hb_mero.items(), key=lambda x: -x[1]):
            f.write(f"  {res[2]} {res[0]}{res[1]}: {cnt}\n")
        f.write("\nMeropenem key residues (Vancomycin-first):\n")
        for res, cnt in sorted(key_hb_mero_van_first.items(), key=lambda x: -x[1]):
            f.write(f"  {res[2]} {res[0]}{res[1]}: {cnt}\n")

    # Save cross-comparison summary
    with open(OUTPUT / "cross_comparison.txt", "w", encoding="utf-8") as f:
        f.write("Cross-comparison summary\n")
        f.write("=" * 60 + "\n")
        f.write(f"Meropenem apo affinity: {e_mero_apo:.3f} kcal/mol\n")
        f.write(f"Vancomycin apo affinity: {e_van_apo:.3f} kcal/mol\n")
        f.write(f"Meropenem in Vancomycin-first complex: {e_van_first_mero:.3f} kcal/mol\n")
        f.write(f"Vancomycin in Meropenem-first complex: {e_mero_first_van:.3f} kcal/mol\n")
        if e_mero_first_van is not None:
            ddG = e_mero_first_van - e_van_first_mero
            f.write(f"DeltaDeltaG (Mero-first - Van-first): {ddG:.3f} kcal/mol\n")
            f.write(f"Vancomycin affinity loss vs apo: {e_mero_first_van - e_van_apo:.3f} kcal/mol\n")
        f.write(f"Meropenem affinity loss vs apo (Van-first): {e_van_first_mero - e_mero_apo:.3f} kcal/mol\n")
        if rmsd_van is not None:
            f.write(f"Vancomycin RMSD (ternary vs apo): {rmsd_van:.3f} A\n")
        if rmsd_mero is not None:
            f.write(f"Meropenem RMSD (Mero-first vs Van-first): {rmsd_mero:.3f} A\n")
        f.write(f"Protein-Vancomycin H-bonds (Mero-first): {hb_total}\n")
        f.write(f"Protein-Meropenem H-bonds (Mero-first): {hb_mero_total}\n")
        f.write(f"Protein-Meropenem H-bonds (Van-first): {hb_mero_van_first_total}\n")


if __name__ == "__main__":
    main()
