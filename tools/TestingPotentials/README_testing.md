# Testing Potentials Suite

Python-based framework for quickly validating MiLaDy ML interatomic potentials by computing key physical properties via LAMMPS.

**Pipeline:** configure → get equilibrium lattice constant + elastic constants → compute defect formation energies → compute migration barriers.

## 1. Global Configuration: `tester.input`

Single file centralizing the setup for all tests:

- **`lammps_executable`** — path to the LAMMPS binary (compiled with `pair_style milady`)
- **`potential_file`** — path to the `.pot` file produced by MiLaDy
- **`lammps_script_block`** — the LAMMPS commands for `pair_style`, `mass`, and `pair_coeff` (multiline, delimited by `"""`)

Example:
```
lammps_executable = /path/to/lmp_mpi
potential_file = /path/to/lammps_bso4_snap1_params.pot
lammps_script_block = """
pair_style      milady
mass       1  55.84500
pair_coeff      * *        pot.fs Mo
"""
```

Every sub-test reads this file to know which potential/executable to use.

## 2. Lattice Parameter & Elastic Constants: `get_a0/`

**Script:** `get_a0/get_a0.py`  
**Local config:** `get_a0/get_a0.input` — provides `lattice_parameter` (initial guess), `structure_type` (bcc/fcc/a15/c15), and `mpi_command`.

### Workflow (`python get_a0.py -m a0`)

1. Reads both `../tester.input` and `get_a0.input`
2. Validates paths and parameters
3. Generates a LAMMPS input file in `get_a0/elastic/LAMMPS/` by injecting the `pair_style` block from `tester.input`
4. Patches `input_elastic.ini_org` with the structure type, lattice guess, LAMMPS executable, and potential path → writes `input_elastic.ini`
5. Calls **Elastic.git** (external elastic-constants code under `elastic/Elastic.git/`) in three phases: `build` → run LAMMPS in each strain directory (bulk, babulk, c11, c44) → `extract`
6. **Output:** optimized lattice parameter and elastic constants (C11, C12, C44, bulk modulus, etc.)

**Cleanup:** `python get_a0.py -m clean`

## 3. Defect Formation Energies

All formation energy scripts share the same `Minimizer` class pattern:

1. Read atomic config (POSCAR) via ASE
2. Scale the cell to the equilibrium lattice parameter
3. Write LAMMPS data file and run LAMMPS minimization via MPI
4. Parse total energy from LAMMPS output
5. Compute formation energy: **E_f = E_defect − (N_defect / N_bulk) × E_bulk**

### 3a. Vacancy Formation: `Ef_small_vac/`

**Script:** `Ef_small_vac/small_vacancy.py`  
**Pre-built configs** in `config/`: `bulk.poscar`, `1vac.poscar`, `2vac.poscar`, `2vac_2nn.poscar`, `2vac_3nn.poscar`

**Run:** `python small_vacancy.py -m formation`  
**Output:** `small_vac.log` with formation energies for mono-vacancy and di-vacancies (1NN, 2NN, 3NN).  
**Cleanup:** `python small_vacancy.py -m clean`

### 3b. SIA Formation: `Ef_sia/`

**Script:** `Ef_sia/sia.py`  
**Pre-built configs** in `config/`: 22 POSCAR files covering 1-SIA through 5-SIA clusters in various orientations (⟨100⟩, ⟨110⟩, ⟨111⟩, octahedral, tetrahedral, mixed).

**Run:** `python sia.py -m formation`  
**Output:** `sia.log` with formation energies for all SIA configurations.  
**Cleanup:** `python sia.py -m clean`

## 4. Migration Barriers

Migration scripts use a `Migration` class that:

1. Minimizes initial and final states with LAMMPS
2. Runs a LAMMPS **NEB** (Nudged Elastic Band) calculation between the two states
3. Extracts the migration energy barrier from the NEB output
4. Uses `CubicSpline` interpolation and `matplotlib` for plotting the energy path

### 4a. Vacancy Migration: `Emig_vac/`

**Script:** `Emig_vac/mig_vacancy.py`  
**Run:** `python mig_vacancy.py -m migration`  
**Cleanup:** `python mig_vacancy.py -m clean`

### 4b. SIA Migration: `Emig_sia/`

**Script:** `Emig_sia/mig_sia.py`  
**Run:** `python mig_sia.py -m migration`  
**Cleanup:** `python mig_sia.py -m clean`

## Summary

| Step | Directory | What it computes |
|------|-----------|-----------------|
| Config | `tester.input` | LAMMPS exe, potential, pair_style block |
| 1 | `get_a0/` | Equilibrium lattice parameter + elastic constants (C_ij) |
| 2 | `Ef_small_vac/` | Vacancy formation energies (1V, 2V-1NN/2NN/3NN) |
| 3 | `Ef_sia/` | SIA formation energies (1–5 SIA, multiple orientations) |
| 4 | `Emig_vac/` | Vacancy migration barrier (NEB) |
| 5 | `Emig_sia/` | SIA migration barrier (NEB) |

## Dependencies

- Python 3 with `numpy`, `ase`, `scipy`, `matplotlib`
- LAMMPS compiled with `pair_style milady` support
- MPI runtime (`mpirun`)
- **Elastic.git** (bundled under `get_a0/elastic/Elastic.git/`)

## Quick Start

1. Edit `tester.input` with your LAMMPS path, potential file, and pair_style block
2. Edit `get_a0/get_a0.input` with lattice parameter guess, structure type, and MPI command
3. Run `cd get_a0 && python get_a0.py -m a0` to get equilibrium lattice constant and elastic constants
4. Update the hardcoded `a0` value in `Ef_small_vac/small_vacancy.py` and `Ef_sia/sia.py`
5. Run `cd Ef_small_vac && python small_vacancy.py -m formation`
6. Run `cd Ef_sia && python sia.py -m formation`
7. Run `cd Emig_vac && python mig_vacancy.py -m migration`
8. Run `cd Emig_sia && python mig_sia.py -m migration`
