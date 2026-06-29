#!/usr/bin/env python3
"""Create cross-comparison figures and update main report with aligned conclusions."""

import csv
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

WORK = Path(r"D:\6UJN_Meropenem_Vancomycin")
REPORT_DIR = WORK / "04_reports"
FIG_DIR = REPORT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

MERO_FIRST_CSV = WORK / "6UJN_Meropenem_Vancomycin.csv"
VAN_FIRST_CSV = WORK / "input" / "6UJN_Vancomycin_Meropenem.csv"
APO_MERO_CSV = WORK / "input" / "6UJN_Meropenem.csv"
APO_VAN_CSV = WORK / "input" / "6UJN_Vancomycin.csv"

MERO_FIRST_COMPLEX = WORK / "md" / "input" / "complex.pdb"
VAN_FIRST_COMPLEX = WORK / "md" / "input_van_first" / "complex.pdb"


def parse_csv_best(path):
    best = None
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                aff = float(row["Binding Affinity"])
            except (KeyError, ValueError):
                continue
            if best is None or aff < best:
                best = aff
    return best


def parse_pdb_atoms(path):
    atoms = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            record = line[:6]
            if record in ("ATOM  ", "HETATM"):
                elem = line[76:78].strip()
                if not elem:
                    elem = line[12:16].strip()[0]
                atoms.append({
                    "atom_name": line[12:16].strip(),
                    "res_name": line[17:20].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "element": elem,
                })
    return atoms


def centroid(atoms):
    n = len(atoms)
    return (sum(a["x"] for a in atoms) / n,
            sum(a["y"] for a in atoms) / n,
            sum(a["z"] for a in atoms) / n)


def is_heavy(atom):
    return atom["element"] != "H"


def main():
    # Parse docking energies
    e_mero_apo = parse_csv_best(APO_MERO_CSV)
    e_van_apo = parse_csv_best(APO_VAN_CSV)
    e_mero_first_van = parse_csv_best(MERO_FIRST_CSV)
    e_van_first_mero = parse_csv_best(VAN_FIRST_CSV)

    # Parse complexes
    mero_first_atoms = parse_pdb_atoms(MERO_FIRST_COMPLEX)
    van_first_atoms = parse_pdb_atoms(VAN_FIRST_COMPLEX)

    def get_components(atoms):
        protein_ca = [a for a in atoms if a["atom_name"] == "CA"]
        van = [a for a in atoms if a["res_name"] == "VAN" and is_heavy(a)]
        mer = [a for a in atoms if a["res_name"] == "MER" and is_heavy(a)]
        return protein_ca, van, mer

    mero_first_ca, mero_first_van, mero_first_mer = get_components(mero_first_atoms)
    van_first_ca, van_first_van, van_first_mer = get_components(van_first_atoms)

    # Compute COM distances
    def com_distance(atoms1, atoms2):
        c1 = centroid(atoms1)
        c2 = centroid(atoms2)
        return math.sqrt((c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2)

    mero_first_dist = com_distance(mero_first_van, mero_first_mer)
    van_first_dist = com_distance(van_first_van, van_first_mer)

    # Figure 1: Affinity comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    categories = ['Meropenem apo', 'Meropenem in\nVancomycin-first', 'Vancomycin apo', 'Vancomycin in\nMeropenem-first']
    values = [e_mero_apo, e_van_first_mero, e_van_apo, e_mero_first_van]
    colors = ['lightcoral', 'coral', 'lightblue', 'steelblue']
    bars = ax.bar(categories, values, color=colors, alpha=0.8, edgecolor='black')
    ax.set_ylabel('Binding affinity (kcal/mol)', fontsize=12)
    ax.set_title('Cross-comparison: sequential vs apo docking affinities', fontsize=14)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.grid(axis='y', alpha=0.3)
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height - 0.3,
                f'{val:.3f}', ha='center', va='top', fontsize=10, fontweight='bold')
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'cross_comparison_affinities.png', dpi=300)
    plt.close(fig)

    # Figure 2: COM distance comparison
    fig, ax = plt.subplots(figsize=(8, 6))
    categories = ['Meropenem-first\n(VAN-MER)', 'Vancomycin-first\n(VAN-MER)']
    distances = [mero_first_dist, van_first_dist]
    colors = ['steelblue', 'coral']
    bars = ax.bar(categories, distances, color=colors, alpha=0.8, edgecolor='black')
    ax.set_ylabel('Ligand centroid-centroid distance (A)', fontsize=12)
    ax.set_title('Cross-comparison: ligand-ligand distances in docked ternary complexes', fontsize=14)
    ax.grid(axis='y', alpha=0.3)
    for bar, val in zip(bars, distances):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height - 2,
                f'{val:.2f} A', ha='center', va='top', fontsize=11, fontweight='bold')
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'cross_comparison_distances.png', dpi=300)
    plt.close(fig)

    # Figure 3: Side-by-side 3D structures
    fig = plt.figure(figsize=(16, 7))

    ax1 = fig.add_subplot(121, projection='3d')
    ax1.scatter([a["x"] for a in mero_first_ca], [a["y"] for a in mero_first_ca],
                [a["z"] for a in mero_first_ca], c='lightgray', s=5, alpha=0.4, label='Protein CA')
    ax1.scatter([a["x"] for a in mero_first_mer], [a["y"] for a in mero_first_mer],
                [a["z"] for a in mero_first_mer], c='red', s=30, alpha=0.8, label='Meropenem')
    ax1.scatter([a["x"] for a in mero_first_van], [a["y"] for a in mero_first_van],
                [a["z"] for a in mero_first_van], c='blue', s=30, alpha=0.8, label='Vancomycin')
    ax1.set_title(f'Meropenem-first (VAN-MER = {mero_first_dist:.2f} A)', fontsize=12)
    ax1.set_xlabel('X (A)')
    ax1.set_ylabel('Y (A)')
    ax1.set_zlabel('Z (A)')
    ax1.legend()

    ax2 = fig.add_subplot(122, projection='3d')
    ax2.scatter([a["x"] for a in van_first_ca], [a["y"] for a in van_first_ca],
                [a["z"] for a in van_first_ca], c='lightgray', s=5, alpha=0.4, label='Protein CA')
    ax2.scatter([a["x"] for a in van_first_mer], [a["y"] for a in van_first_mer],
                [a["z"] for a in van_first_mer], c='red', s=30, alpha=0.8, label='Meropenem')
    ax2.scatter([a["x"] for a in van_first_van], [a["y"] for a in van_first_van],
                [a["z"] for a in van_first_van], c='blue', s=30, alpha=0.8, label='Vancomycin')
    ax2.set_title(f'Vancomycin-first (VAN-MER = {van_first_dist:.2f} A)', fontsize=12)
    ax2.set_xlabel('X (A)')
    ax2.set_ylabel('Y (A)')
    ax2.set_zlabel('Z (A)')
    ax2.legend()

    fig.tight_layout()
    fig.savefig(FIG_DIR / 'cross_comparison_3d_structures.png', dpi=300)
    plt.close(fig)

    # Figure 4: Energy + distance summary panel
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: affinity comparison
    ax = axes[0]
    categories = ['Meropenem apo', 'Mero in Van-first', 'Vancomycin apo', 'Van in Mero-first']
    values = [e_mero_apo, e_van_first_mero, e_van_apo, e_mero_first_van]
    colors = ['lightcoral', 'coral', 'lightblue', 'steelblue']
    bars = ax.bar(categories, values, color=colors, alpha=0.8, edgecolor='black')
    ax.set_ylabel('Binding affinity (kcal/mol)', fontsize=11)
    ax.set_title('Docking affinity comparison', fontsize=12)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.grid(axis='y', alpha=0.3)
    for bar, val in zip(bars, values):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height - 0.3,
                f'{val:.3f}', ha='center', va='top', fontsize=9, fontweight='bold')

    # Right: COM distance comparison
    ax = axes[1]
    categories = ['Mero-first', 'Van-first']
    distances = [mero_first_dist, van_first_dist]
    colors = ['steelblue', 'coral']
    bars = ax.bar(categories, distances, color=colors, alpha=0.8, edgecolor='black')
    ax.set_ylabel('VAN-MER centroid distance (A)', fontsize=11)
    ax.set_title('Docked ternary complex geometry', fontsize=12)
    ax.grid(axis='y', alpha=0.3)
    for bar, val in zip(bars, distances):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height - 2,
                f'{val:.2f} A', ha='center', va='top', fontsize=10, fontweight='bold')

    fig.suptitle('Cross-comparison summary: docking energy vs ternary geometry', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(FIG_DIR / 'cross_comparison_summary.png', dpi=300)
    plt.close(fig)

    print(f"Saved cross-comparison figures to {FIG_DIR}")

    # Save cross-comparison data
    summary_path = REPORT_DIR / "cross_comparison_docking_data.txt"
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("Cross-comparison: docking energies and geometries\n")
        f.write("=" * 60 + "\n\n")
        f.write("Docking energies (kcal/mol):\n")
        f.write(f"  Meropenem apo:             {e_mero_apo:.3f}\n")
        f.write(f"  Meropenem in Vancomycin-first: {e_van_first_mero:.3f}\n")
        f.write(f"  Vancomycin apo:            {e_van_apo:.3f}\n")
        f.write(f"  Vancomycin in Meropenem-first: {e_mero_first_van:.3f}\n\n")
        f.write("Docked ternary complex geometries:\n")
        f.write(f"  Meropenem-first VAN-MER distance: {mero_first_dist:.2f} A\n")
        f.write(f"  Vancomycin-first VAN-MER distance: {van_first_dist:.2f} A\n\n")
        f.write("MM-GBSA binding free energies (igb=2, saltcon=0.15 M):\n")
        f.write("  Meropenem-first (Vancomycin as ligand):  -52.8460 +/- 4.2025 kcal/mol\n")
        f.write("  Vancomycin-first (Meropenem as ligand):  -23.0253 +/- 2.4839 kcal/mol\n")
        f.write("  Delta Delta G (Mero-first - Van-first):  -29.82 kcal/mol\n\n")
        f.write("Interpretation:\n")
        f.write("  - Meropenem-first: Vancomycin binds close to Meropenem, forming a crowded ternary complex.\n")
        f.write("  - Vancomycin-first: Meropenem binds far from Vancomycin, no ternary crowding; affinity is essentially apo-like.\n")
        f.write("  - Docking affinities alone suggest little difference, but geometry + MM-GBSA together reveal order-dependent paving.\n")

    print(f"Saved cross-comparison data: {summary_path}")


if __name__ == "__main__":
    main()
