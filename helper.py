import profile
import time
import re
from pathlib import Path
from ase.io import write, read
from ase import Atoms   
from ase.calculators.espresso import Espresso, EspressoProfile
import sys
import subprocess
from threading import Thread
from collections import Counter

## Pathways
base_path = Path("/Users/famkepotze/Desktop/3S/DFT_folder")

script_path = base_path / "NV_center_QE/scripts"
results_path = base_path / "results"

structures_path = results_path / "structures"
calculations_path = results_path / "calculations"


qe_path   = Path("/Users/famkepotze/Desktop/3S/NV_center_3S/QE_extra_files/qe-7.5")

pw_path   = qe_path / "bin/pw.x"
pseudo_dir = qe_path / "pseudo"
dos_path  = qe_path / "bin/dos.x"
projwfc_path = qe_path / "bin/projwfc.x"



# ─────────────────────────────────────────────────────────────
# Material-specific settings
# ─────────────────────────────────────────────────────────────

## variables
# ─────────────────────────────────────────────────────────────
supercell_sizes = {
    1: (1, 1, 1),
    2: (2, 2, 2),
    3: (3, 3, 3)
}

supercell_labels = {
    1: "1x1x1",    
    2: "2x2x2",
    3: "3x3x3"
}

materials = [
    'diamond', 
    'NV_center',
    ]


# k-points scaled down with supercell size to keep reciprocal-space
# sampling roughly constant: 6x6x6 for 1x1x1 → 3x3x3 for 2x2x2 → 2x2x2 for 3x3x3
k_points_per_supercell = {
    1: (6, 6, 6),
    2: (3, 3, 3),
    3: (2, 2, 2),
}

## QE Profile and Pseudopotentials
profile = EspressoProfile(
    command=str(pw_path),
    pseudo_dir=str(pseudo_dir)
)

pseudopotentials = {
    'C': 'C.pbe-n-rrkjus_psl.1.0.0.UPF',
    'N': 'N.pbe-n-rrkjus_psl.1.0.0.UPF'
}


def get_material_settings(material):
    """Return system card settings for a given material."""
    if material == "diamond":
        return {
            "nspin": 1,
            "tot_charge": 0,
            # no tot_magnetization needed for nspin=1
            "occupations": "fixed",
        }
    elif material == "NV_center":
        return {
            "nspin": 2,
            "tot_magnetization": 2,
            "tot_charge": -1,
            "occupations": "smearing",
            "smearing": "mv",
            "degauss": 0.01,
        }
    else:
        raise ValueError(f"Unknown material: {material!r}. "
                         "Add its settings to get_system_settings().")




from pathlib import Path
from ase.io import read
from ase import Atoms

from pathlib import Path
from ase.io import read
from ase import Atoms
import subprocess


class QEWorkflowEngine:

    def __init__(self):

        self.labels = {
            1: "1x1x1",
            2: "2x2x2",
            3: "3x3x3",
        }

        self.k_scf = {
            1: (7, 7, 7),
            2: (2, 2, 2),
            3: (1, 1, 1),
        }

        self.k_nscf = {
            1: (12, 12, 12),
            2: (6, 6, 6),
            3: (4, 4, 4),
        }

        self.base = base_path
        self.relaxed = structures_path / "relaxed_structures"
        self.out = self.base / calculations_path
        self.out.mkdir(parents=True, exist_ok=True)

    # ======================================================
    # STRUCTURE
    # ======================================================

    def load(self, material, size):

        label = self.labels[size]
        path = self.relaxed / f"{material}_{label}.xyz"

        raw = read(path)

        return Atoms(
            symbols=raw.get_chemical_symbols(),
            positions=raw.get_positions(),
            cell=raw.get_cell(),
            pbc=True
        )

    # ======================================================
    # SETTINGS
    # ======================================================

    def settings(self, material):

        if material == "diamond":
            return dict(
                nspin=1,
                occupations="fixed",
                charge=0
            )

        elif material == "NV_center":
            return dict(
                nspin=2,
                occupations="smearing",
                charge=-1,
                degauss=0.0005
            )

    # ======================================================
    # BUILDERS
    # ======================================================

    def control(self, prefix, calc):

        return f"""&control
    calculation = '{calc}',
    prefix = '{prefix}',
    outdir = './out',
    pseudo_dir = '{str(pseudo_dir)}',
/
"""

    def electrons(self):

        return """&electrons
    conv_thr = 1.0d-8,
    mixing_beta = 0.3,
/
"""

    def system(self, structure, material):

        s = self.settings(material)

        nat = len(structure)
        ntyp = len(set(structure.get_chemical_symbols()))

        block = f"""&system
    ibrav = 0,
    nat = {nat},
    ntyp = {ntyp},
    ecutwfc = 60,
    ecutrho = 480,
    nspin = {s['nspin']},
"""

        if material == "NV_center":
            block += f"""    tot_charge = {s['charge']},
    occupations = '{s['occupations']}',
    smearing = 'gaussian',
    degauss = {s['degauss']},
    starting_magnetization(1) = 0.0,
    starting_magnetization(2) = 0.5,
"""
        else:
            block += f"""    occupations = '{s['occupations']}',
"""

        block += "/\n"
        return block

    def atoms(self, structure):

        return "\n".join(
            f"{a.symbol} {a.position[0]:.6f} {a.position[1]:.6f} {a.position[2]:.6f}"
            for a in structure
        )

    def cell(self, structure):

        return "\n".join(
            f"{v[0]:.6f} {v[1]:.6f} {v[2]:.6f}"
            for v in structure.cell.array
        )

    def species(self, structure):

        syms = set(structure.get_chemical_symbols())
        lines = []

        if "C" in syms:
            lines.append("C 12.011 C.pbe-n-rrkjus_psl.1.0.0.UPF")

        if "N" in syms:
            lines.append("N 14.007 N.pbe-n-rrkjus_psl.1.0.0.UPF")

        return "\n".join(lines)

    def kpoints(self, size, mode="scf"):

        k = self.k_scf[size] if mode == "scf" else self.k_nscf[size]

        return f"""K_POINTS automatic
{k[0]} {k[1]} {k[2]} 0 0 0
"""

    # ======================================================
    # WRITER
    # ======================================================

    def write_scf(self, material, size):

        structure = self.load(material, size)
        label = self.labels[size]
        prefix = f"{material}_{label}"

        folder = self.out / material / label
        folder.mkdir(parents=True, exist_ok=True)

        text = (
            self.control(prefix, "scf")
            + self.system(structure, material)
            + self.electrons()
            + "\nCELL_PARAMETERS angstrom\n"
            + self.cell(structure)
            + "\n\nATOMIC_SPECIES\n"
            + self.species(structure)
            + "\n\nATOMIC_POSITIONS angstrom\n"
            + self.atoms(structure)
            + "\n\n"
            + self.kpoints(size, "scf")
        )

        (folder / "scf.in").write_text(text)

    def write_nscf(self, material, size):

        structure = self.load(material, size)
        label = self.labels[size]
        prefix = f"{material}_{label}"

        folder = self.out / material / label

        text = (
            self.control(prefix, "nscf")
            + self.system(structure, material)
            + self.electrons()
            + "\nCELL_PARAMETERS angstrom\n"
            + self.cell(structure)
            + "\n\nATOMIC_SPECIES\n"
            + self.species(structure)
            + "\n\nATOMIC_POSITIONS angstrom\n"
            + self.atoms(structure)
            + "\n\n"
            + self.kpoints(size, "nscf")
        )

        (folder / "nscf.in").write_text(text)

    def write_dos(self, material, size):

        label = self.labels[size]
        prefix = f"{material}_{label}"

        folder = self.out / material / label

        text = f"""&dos
    prefix = '{prefix}',
    outdir = './out',
    fildos = '{prefix}.dos',
    Emin = -20.0,
    Emax = 40.0,
    DeltaE = 0.01,
/
"""

        (folder / "dos.in").write_text(text)

    def write_projwfc(self, material, size):

        label = self.labels[size]
        prefix = f"{material}_{label}"

        folder = self.out / material / label

        text = f"""&projwfc
    prefix = '{prefix}',
    outdir = './out',
    filpdos = '{prefix}.pdos',
/
"""

        (folder / "projwfc.in").write_text(text)

# ============================================================
# IMPORTS
# ============================================================

import time
import subprocess
import numpy as np

from pathlib import Path
from threading import Thread

from ase.io import read, write
from ase.calculators.espresso import Espresso


 
# ============================================================
# CLASS: StructureRelaxation
# ============================================================
 
class StructureRelaxation:
 
    def __init__(self):
 
        self.labels = supercell_labels         # {1: "1x1x1", ...}
        self.kpts   = k_points_per_supercell   # {1: (6,6,6), ...}
 
    # ========================================================
    # PATHS
    # ========================================================
 
    def get_input_structure(self, material, size):
 
        label = self.labels[size]
 
        return (
            structures_path /
            "raw_structures" /
            f"{material}_{label}.xyz"
        )
 
    def get_output_structure(self, material, size):
 
        label  = self.labels[size]
        output = (
            structures_path /
            "relaxed_structures" /
            f"{material}_{label}.xyz"
        )
        output.parent.mkdir(parents=True, exist_ok=True)
 
        return output
 
    def get_calculation_directory(self, material, size):
 
        label  = self.labels[size]
        folder = calculations_path / material / label
        folder.mkdir(parents=True, exist_ok=True)
 
        return folder
 
    # ========================================================
    # LOAD STRUCTURE
    # ========================================================
 
    def load_structure(self, material, size):
 
        return read(self.get_input_structure(material, size))
 
    # ========================================================
    # CREATE CALCULATOR
    # ========================================================
 
    def create_calculator(self, material, size, calculation_directory):
 
        calc = Espresso(
 
            profile=profile,
 
            pseudopotentials=pseudopotentials,
 
            input_data={
 
                "control": {
                    "calculation": "relax",
                    "pseudo_dir": str(pseudo_dir),
                    "outdir": "./out",
                    "tprnfor": True,
                    "tstress": True,
                },
 
                "system": {
                    "ecutwfc": 60,
                    "ecutrho": 480,
                    **get_material_settings(material),
                },
 
                "electrons": {
                    "conv_thr": 1e-8,
                    "mixing_beta": 0.3,
                },
 
                "ions": {
                    "ion_dynamics": "bfgs",
                },
            },
 
            kpts=self.kpts[size],
 
            directory=str(calculation_directory),
        )
 
        return calc
 
    # ========================================================
    # RUN SINGLE RELAXATION
    # ========================================================
 
    def run_relaxation(self, material, size):
 
        output_structure      = self.get_output_structure(material, size)
        calculation_directory = self.get_calculation_directory(material, size)
        atoms                 = self.load_structure(material, size)
        pwo_path              = calculation_directory / "espresso.pwo"
 
        monitor = ProgressUpdate(pwo_path, mode="relax")
        t       = Thread(target=monitor.run)
        t.start()
 
        atoms.calc = self.create_calculator(material, size, calculation_directory)
 
        start_time = time.time()
 
        try:
 
            energy    = atoms.get_potential_energy()
            forces    = atoms.get_forces()
            max_force = np.max(np.linalg.norm(forces, axis=1))
            runtime   = time.time() - start_time
 
            t.join()
 
            print(f"\n  {'='*44}")
            print(f"  RELAXATION FINISHED")
            print(f"  Final energy : {energy:.6f} eV")
            print(f"  Max force    : {max_force:.6f} eV/Å")
            print(f"  Runtime      : {runtime:.2f} s")
            print(f"  Saved to     : {output_structure}")
            print(f"  {'='*44}")
 
            write(output_structure, atoms)
 
        except Exception as e:
 
            print(f"\n  {'='*44}")
            print(f"  CALCULATION FAILED")
            print(f"  Error: {e}")
            print(f"  {'='*44}")
 
            if pwo_path.exists():
 
                print("\n  --- LAST 200 LINES OF QE OUTPUT ---\n")
 
                result = subprocess.run(
                    ["tail", "-n", "200", str(pwo_path)],
                    capture_output=True,
                    text=True,
                )
 
                print(result.stdout)




# ─────────────────────────────────────────────────────────────
# Progress monitor
# ─────────────────────────────────────────────────────────────
class ProgressUpdate:
    def __init__(self, pwo_path, mode="scf", interval=5):
        self.pwo_path   = Path(pwo_path)
        self.mode       = mode
        self.interval   = interval
        self.last_scf_iter = -1
        self.last_energy   = None
        self.ionic_step    = 0

    def parse(self, text):
        lines = text.splitlines()
        for line in reversed(lines):
            if "iteration #" in line:
                try:
                    it = int(line.split()[2])
                    if it != self.last_scf_iter:
                        self.last_scf_iter = it
                        print(f"[LOGGER] Iteration: {it}")
                    break
                except:
                    pass
        for line in reversed(lines):
            if "!    total energy" in line:
                try:
                    energy = float(line.split()[-2])
                    if self.last_energy is not None:
                        dE = energy - self.last_energy
                        print(f"[LOGGER] Energy: {energy:.6f}  ΔE={dE:.2e}")
                    else:
                        print(f"[LOGGER] Energy: {energy:.6f}")
                    self.last_energy = energy
                    break
                except:
                    pass
        if self.mode == "relax":
            ionic_steps = len(re.findall(r"Entering BFGS Geometry Optimization", text))
            if ionic_steps != self.ionic_step:
                self.ionic_step = ionic_steps
                print(f"[LOGGER] Ionic step: {self.ionic_step}")

    def run(self):
        print("[LOGGER] QE Progress Monitor started...\n")
        while True:
            if self.pwo_path.exists():
                with open(self.pwo_path, "r") as f:
                    text = f.read()
                self.parse(text)
                if "convergence has been achieved" in text:
                    print("\n[LOGGER] SCF Converged!")
                    break
                if "End of BFGS Geometry Optimization" in text:
                    print("\n[LOGGER] Relaxation Finished!")
                    break
            time.sleep(self.interval)


# ─────────────────────────────────────────────────────────────
