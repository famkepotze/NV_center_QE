
from helper import (
    QEWorkflowEngine,
    calculations_path,
)

engine = QEWorkflowEngine()

materials = [
    "diamond",
    "NV_center"
]

sizes = [1, 2, 3]

skip = {
    ("diamond", 2),
    ("diamond", 3),
}


print("\n==================================================")
print("GENERATING QE INPUT FILES")
print("==================================================")


for material in materials:

    for size in sizes:

        label = engine.labels[size]

        if (material, size) in skip:
            print(f"\n  {label}  ← SKIPPED")
            continue

        folder = (
            calculations_path /
            material /
            label
        )

        # --------------------------------------------------
        # WRITE FILES
        # --------------------------------------------------

        engine.write_scf(material, size)

        engine.write_nscf(material, size)

        engine.write_dos(material, size)

        engine.write_projwfc(material, size)

        # --------------------------------------------------
        # PRINT INFO
        # --------------------------------------------------

        print("\n--------------------------------------------------")

        print(
            f"Generated QE inputs for:\n"
            f"{material} | {label}"
        )

        print(f"\nFolder:\n{folder}")

        print("\nFiles:")

        print(folder / "scf.in")

        print(folder / "nscf.in")

        print(folder / "dos.in")

        print(folder / "projwfc.in")


print("\n==================================================")
print("ALL INPUT FILES GENERATED")
print("==================================================")
