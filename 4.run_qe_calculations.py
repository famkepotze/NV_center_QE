import subprocess
import time
from datetime import datetime

from helper import (
    calculations_path,
    pw_path,
    dos_path,
    projwfc_path,
)

materials = [
    #"diamond", 
    "NV_center"
    ]

sizes = [
    #"1x1x1",
    #"2x2x2",
    "3x3x3",
]

skip = {
    ("diamond", "2x2x2"),
    ("diamond", "3x3x3"),
}


# ======================================================
# HELPER
# ======================================================

def timestamp():

    return datetime.now().strftime(
        "%Y-%m-%d %H:%M:%S"
    )

def runtime_string(seconds):

    hours = int(seconds // 3600)

    minutes = int((seconds % 3600) // 60)

    seconds = int(seconds % 60)

    return (
        f"{hours:02d}h "
        f"{minutes:02d}m "
        f"{seconds:02d}s"
    )


def run(cmd, cwd, step_name):

    print("\n--------------------------------------------------")
    print(f"{step_name}")
    print("--------------------------------------------------")

    print(f"Start : {timestamp()}")
    print(f"Folder: {cwd}")
    print(f"Command:\n{cmd}")

    start = time.time()

    subprocess.run(
        cmd,
        shell=True,
        cwd=cwd,
        check=True
    )

    end = time.time()

    runtime = end - start

    print(f"Finish: {timestamp()}")
    print(
        f"Runtime: "
        f"{runtime_string(runtime)}"
    )


# ======================================================
# MAIN LOOP
# ======================================================

workflow_start = time.time()

print("\n==================================================")
print("STARTING QUANTUM ESPRESSO WORKFLOW")
print("==================================================")

print(f"Global start time: {timestamp()}")


for material in materials:

    for size in sizes:

        folder = (
            calculations_path /
            material /
            size
        )

        label = material, size

        if (material, size) in skip:
            print(f"\n  {label}  ← SKIPPED")
            continue


        print("\n==================================================")
        print(f"{material}  |  {size}")
        print("==================================================")

        # ---------------- SCF ----------------

        run(
            f"{pw_path} < scf.in > scf.out",
            folder,
            "SCF"
        )

        # ---------------- NSCF ----------------

        run(
            f"{pw_path} < nscf.in > nscf.out",
            folder,
            "NSCF"
        )

        # ---------------- DOS ----------------

        run(
            f"{dos_path} < dos.in > dos.out",
            folder,
            "DOS"
        )

        # ---------------- PDOS ----------------

        run(
            (
                f"{projwfc_path} "
                f"< projwfc.in > projwfc.out"
            ),
            folder,
            "PROJWFC"
        )


workflow_end = time.time()

print("\n==================================================")
print("ALL CALCULATIONS FINISHED")
print("==================================================")

print(f"Global finish time: {timestamp()}")

print(
    "Total workflow runtime: "
    f"{runtime_string(workflow_end - workflow_start)}"
)
