# NV Center — DFT Workflow

A computational workflow for calculating the electronic structure and density of states (DOS/PDOS) of diamond and nitrogen-vacancy (NV⁻) center supercells using Quantum ESPRESSO (QE) and ASE.

---

## Overview

This project automates the full DFT pipeline from structure building to DOS plotting across three supercell sizes (1×1×1, 2×2×2, 3×3×3) for two materials:

- **Diamond** — pristine carbon supercell, nspin=1, fixed occupations
- **NV⁻ center** — one C replaced by N, one C vacancy, charge=-1, nspin=2

The pipeline runs: structure building → geometry relaxation → SCF → NSCF → DOS → PDOS → plotting.

---

## Project Structure

```
NV_center_QE/
│
├── scripts/
│   ├── helper.py                  # All classes, paths, and settings
│   ├── 1_build_structure.ipynb    # Build raw supercell structures
│   ├── 2_structure_relaxation.py  # Geometry relaxation via ASE + QE
│   ├── 3_write_qe_input.py        # Generate SCF/NSCF/DOS/PDOS input files
│   ├── 4_run_qe_calculations.py   # Run QE executables in sequence
│   └── 5_plot_results.ipynb       # Plot DOS and PDOS results
│
└── results/
    ├── structures/
    │   ├── raw_structures/        # Input .xyz files from step 1
    │   └── relaxed_structures/    # Relaxed .xyz files from step 2
    └── calculations/
        ├── diamond/
        │   └── 1x1x1/             # scf.in, nscf.in, dos.in, projwfc.in, out/
        └── NV_center/
            ├── 1x1x1/
            ├── 2x2x2/
            └── 3x3x3/
```

---

## Requirements

- Python 3.8+
- [ASE](https://wiki.fysik.dtu.dk/ase/) — `pip install ase`
- [Quantum ESPRESSO 7.5](https://www.quantum-espresso.org/) — compiled and accessible at `qe-7.5/bin/`
- Pseudopotentials (PBE, RRKJUS):
  - `C.pbe-n-rrkjus_psl.1.0.0.UPF`
  - `N.pbe-n-rrkjus_psl.1.0.0.UPF`

Set your base path in `helper.py`:
```python
base_path = Path("/path/to/your/DFT_folder")
```

---

## Workflow

### Step 1 — Build structures (`1_build_structure.ipynb`)
Constructs diamond and NV center supercells for sizes 1×1×1, 2×2×2, and 3×3×3. For the NV center, one carbon atom is substituted with nitrogen and the adjacent carbon is removed. Structures are saved as `.xyz` files in `results/structures/raw_structures/`.

### Step 2 — Geometry relaxation (`2_structure_relaxation.py`)
Runs a `relax` calculation via ASE's `Espresso` calculator using BFGS ion dynamics. A `ProgressUpdate` thread monitors the QE output live, printing SCF iterations, energy, and ionic steps. Relaxed structures are saved to `results/structures/relaxed_structures/`.

```bash
python 2_structure_relaxation.py
```

### Step 3 — Write QE input files (`3_write_qe_input.py`)
Uses `QEWorkflowEngine` to generate four input files per material/size combination: `scf.in`, `nscf.in`, `dos.in`, `projwfc.in`. The number of bands (`nbnd`) is calculated automatically from the valence electron count of the structure.

```bash
python 3_write_qe_input.py
```

### Step 4 — Run QE calculations (`4_run_qe_calculations.py`)
Runs the full SCF → NSCF → DOS → PDOS sequence by calling QE executables via `subprocess`. Outputs are written to `out/` inside each calculation folder.

```bash
python 4_run_qe_calculations.py
```

### Step 5 — Plot results (`5_plot_results.ipynb`)
Reads `.dos` and `.pdos` output files and plots the total and projected density of states.

---

## Material Settings

| Parameter | Diamond | NV⁻ Center |
|---|---|---|
| `nspin` | 1 | 2 |
| `tot_charge` | 0 | -1 |
| `occupations` (SCF) | `fixed` | `smearing` (gaussian) |
| `occupations` (NSCF) | `fixed` | `tetrahedra` |
| `degauss` | — | 0.0005 Ry |
| `starting_magnetization(C)` | — | 0.0 |
| `starting_magnetization(N)` | — | 0.5 |

---

## Supercell Settings

| Size | Atoms (diamond) | Atoms (NV⁻) | k-points SCF | k-points NSCF |
|---|---|---|---|---|
| 1×1×1 | 8 | 7 | 7×7×7 | 12×12×12 |
| 2×2×2 | 64 | 63 | 2×2×2 | 6×6×6 |
| 3×3×3 | 216 | 215 | 1×1×1 | 4×4×4 |

`nbnd` is computed per structure as `max(n_occ + 20, n_occ * 1.2)` where `n_occ` is derived from the valence electron count.

---

## Key Classes (`helper.py`)

**`QEWorkflowEngine`** — generates all QE input files. Methods:
- `write_scf(material, size)` — SCF input
- `write_nscf(material, size)` — NSCF input (tetrahedra for NV center)
- `write_dos(material, size)` — DOS input
- `write_projwfc(material, size)` — PDOS input

**`StructureRelaxation`** — runs geometry relaxation via ASE. Methods:
- `run_relaxation(material, size)` — full relax calculation with live progress output

**`ProgressUpdate`** — background thread that tails the QE `.pwo` output and prints SCF iteration, energy convergence, and ionic steps in real time.

---

## Notes

- Diamond 2×2×2 and 3×3×3 are skipped by default in steps 3 and 4 (only 1×1×1 is needed as reference). Toggle via the `skip` set in each script.
- Always clear the `out/` directory between runs with different prefixes to avoid QE reading stale charge densities.
- The PBE functional underestimates the diamond band gap (~4.1 eV vs experimental 5.47 eV). This is expected behavior for GGA.