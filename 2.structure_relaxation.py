import numpy as np
import subprocess
import time
from threading import Thread
from ase.io import write, read
from ase.calculators.espresso import Espresso

from helper import (
    StructureRelaxation,
    calculations_path,
    structures_path,
    supercell_labels,
    k_points_per_supercell,
)

# ============================================================
# CONFIGURATION
# ============================================================

relaxer = StructureRelaxation()

materials = [
    "diamond",
    "NV_center",
]

sizes = [1, 2, 3]

# Combinations to skip — add/remove as needed
skip = {
    ("diamond", 2),
    ("diamond", 3),
}

labels = supercell_labels          # {1: "1x1x1", 2: "2x2x2", 3: "3x3x3"}
kpts   = k_points_per_supercell    # {1: (6,6,6),  2: (3,3,3),  3: (2,2,2)}

# ============================================================
# RUN RELAXATIONS
# ============================================================

print("\n==================================================")
print("STARTING QE STRUCTURE RELAXATIONS")
print("==================================================")
print(f"Materials : {materials}")
print(f"Sizes     : {[labels[s] for s in sizes]}")

total_start = time.time()

for material in materials:

    print(f"\n{'='*50}")
    print(f"MATERIAL: {material}")
    print(f"{'='*50}")

    for size in sizes:

        label  = labels[size]

        if (material, size) in skip:
            print(f"\n  {label}  ← SKIPPED")
            continue

        folder = calculations_path / material / label

        input_structure  = structures_path / "raw_structures"     / f"{material}_{label}.xyz"
        output_structure = structures_path / "relaxed_structures"  / f"{material}_{label}.xyz"
        
        if output_structure.exists():
            print(f"\n  {label}  ← ALREADY DONE ({output_structure})")
            continue
        
        # --------------------------------------------------
        # PRINT INFO
        # --------------------------------------------------

        print(f"\n  {label}")
        print(f"  k-points        : {kpts[size]}")
        print(f"  Input structure : {input_structure}")
        print(f"  Calc folder     : {folder}")
        print(f"  Output          : {output_structure}")

        # --------------------------------------------------
        # RUN
        # --------------------------------------------------

        relaxer.run_relaxation(material, size)


total_runtime = time.time() - total_start

print("\n==================================================")
print("ALL RELAXATIONS FINISHED")
print(f"Total runtime : {total_runtime / 60:.2f} minutes")
print("==================================================")