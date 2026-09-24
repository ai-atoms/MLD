# MILADY

[![docs build](https://img.shields.io/github/actions/workflow/status/ai-atoms/milady-docs/main.yml?branch=main&label=docs%20build&style=for-the-badge)](https://ai-atoms.github.io/milady-docs/)
[![licence](https://img.shields.io/badge/licence-ASL-blue?style=for-the-badge)](LICENSE.md)
[![PyPI](https://img.shields.io/badge/PyPI-milady--mlip-orange?style=for-the-badge)](https://pypi.org/project/milady-mlip/)

**MILADY** (*Machine Learning Dynamics*) builds machine-learning interatomic
potentials and regression models from *ab initio* data. It aims to make
atomistic simulations more accurate and predictive while keeping the
computational cost reasonable. Potentials trained with MILADY run in
[LAMMPS](https://www.lammps.org) for large-scale molecular dynamics.

📖 **Full documentation: <https://ai-atoms.github.io/milady-docs/>**

---

## Contents

- [Features](#features)
- [Installation](#installation)
- [Running MILADY](#running-milady)
- [Using MILADY potentials in LAMMPS](#using-milady-potentials-in-lammps)
- [Examples](#examples)
- [Tests](#tests)
- [Repository layout](#repository-layout)
- [Publications](#publications)
- [Contributors](#contributors)
- [Licence](#licence)

## Features

**Main functionalities**

- Train machine-learning force fields on energies, forces and stresses for
  molecular dynamics. 
- Represent atomic structures in the feature space of atomic descriptors, for
  structural analysis and visualisation.
- An optional physics-based short-range **ZBL** repulsive
  core and an explicit two-body **k2b** channel handle close interatomic
  distances.
- Analyse, sparsify and optimise datasets, for example by choosing reference
  points with MCD/Mahalanobis or CUR selection.
- Read training databases in several formats:
  - a directory of files in the native `.poscar` format (VASP-compatible, with
    energy, forces, stress and spin tags);
  - a single **extended XYZ** file (`.xyz` / `.extxyz`), the format used by
    [ASE](https://wiki.fysik.dtu.dk/ase/);
  - a single **MPtrj-style JSON** file (`.json`).

  ASE `.traj` files can be exported to `.extxyz` or `.poscar` with ASE.

**Descriptors.** MILADY has a large library of atomic descriptors. For the
details, see the
[descriptors page](https://ai-atoms.github.io/milady-docs/contents/ml/descriptors.html).

**Models.** MILADY offers a range of regression models that are linear in the
parameters, including kernel models. For the details, see the
[models page](https://ai-atoms.github.io/milady-docs/contents/ml/model.html)
and the [kernels page](https://ai-atoms.github.io/milady-docs/contents/ml/kernels.html).

**Built for HPC**

- Parallelised with MPI, PBLAS and ScaLapack.
- Scales to large training sets. A streaming mode avoids storing the full
  design matrix.

## Installation

### From PyPI

```bash
pip install milady-mlip
```

### From source (CMake)

CMake is the only supported and tested way to build MILADY.

**Prerequisites**

| Dependency | Tested versions |
|------------|-----------------|
| Fortran compiler | `ifort` > 2020, `ifx` > 2023, `gfortran` > 12 |
| LAPACK / ScaLapack | Intel MKL > 19 |
| MPI | OpenMPI > 4 or Intel MPI > 2020 |
| CMake | > 3.24 |
| HDF5, Fortran interface (*optional*) | 1.14.3 is bundled |

The detailed build steps (environment variables, the Intel/GNU/mixed build
modes, optional HDF5 support) are on the
[local build page](https://ai-atoms.github.io/milady-docs/contents/install/installation.html#local-build).
Once built, the executable is `bin/milady_main.exe`.

The documentation also has instructions for
[HPC clusters (Irene, Marconi)](https://ai-atoms.github.io/milady-docs/contents/install/install_irene.html)
and the [former Makefile build](https://ai-atoms.github.io/milady-docs/contents/install/install_legacy.html).

## Running MILADY

A run directory needs four inputs:

| File | Purpose |
|------|---------|
| `name.in` | One line with the simulation prefix `PREFIXSIM` |
| `PREFIXSIM.ml` | Main input file (`&input_ml` namelist): model, descriptors, database, weight optimisation |
| `db_model.in` | Describes how the database is used: classes, train/test split, energy/force/stress tags and weights |
| The dataset | A `DB/` directory of `.poscar` files (the default), or a single extended XYZ (`.xyz` / `.extxyz`) or JSON (`.json`) file, chosen with `db_path` |

See the [database page](https://ai-atoms.github.io/milady-docs/contents/ml/database.html)
for the file formats and the `db_model.in` syntax.

```bash
mpirun -np 4 milady_main.exe > out.run
```

Here is a minimal `PREFIXSIM.ml` for a linear fit of bcc Fe with the bispectrum
SO(4) descriptor:

```fortran
&input_ml
  ml_type           = 0        ! regression with basis functions
  mld_order         = 1        ! 1 = LML, 2 = quadratic, 3 = polynomial chaos
  desc_forces       = .true.   ! fit forces
  weighted          = .false.  ! .true. for multi-component systems
  chemical_elements = "Fe"
  r_cut             = 4.7d0    ! cutoff radius (Å)
  descriptor_type   = 9        ! 9 = bispectrum SO(4), 4 = AFS, 300 = ACE/kACE
  j_max             = 4.0
&end
```

`ml_type` selects the task: `0` is regression with basis functions, `1` is
kernel ridge regression, `-1` computes descriptors only, and `-2` analyses the
data and chooses the kernel. See the
[input reference](https://ai-atoms.github.io/milady-docs/contents/ml/input.html)
for every keyword.

### From Python

The PyPI package includes a small wrapper:

```python
from milady import milady_run

res = milady_run(
    input="./input.ml",
    db_model="./db_model.in",
    DB="./DB",
    nprocs=2,
    workdir="./",
    output="out.run",   # None: no extra log (milady.out is always written)
    mpi_args=[],
    check=True,         # raise MiladyError on a non-zero exit code
)
```

## Using MILADY potentials in LAMMPS

Build LAMMPS with the `ml-milady` library from
[milady_lammps](https://github.com/mcmarinica/milady_lammps), following the
[LAMMPS installation guide](https://ai-atoms.github.io/milady-docs/contents/install/install_milady-lammps.html).
Then, in your LAMMPS input:

```lammps
pair_style milady
pair_coeff * * Fe_LML.pot Fe
# multi-element, e.g. a TaTiVW high-entropy alloy:
# pair_coeff * * hea.pot Ta Ti V W

# optional: per-atom descriptors
compute D all milady/atom
dump 1 all custom 10 dump.milady.* id type c_D[*]
```

Ready-to-use W and Fe potentials (LML, QNML, KNML) are listed on the
[potentials page](https://ai-atoms.github.io/milady-docs/contents/pubs/potentials.html).
The [`tools/TestingPotentials`](tools/TestingPotentials) suite checks a
potential in LAMMPS by computing its lattice constant, elastic constants,
defect formation energies and migration barriers.

## Examples

Each folder in [`examples/`](examples) is a complete run directory. The
matching documentation page explains it step by step.

| # | Folder | What it shows |
|---|--------|---------------|
| 1 | [`compute_desc_only`](examples/compute_desc_only) | Compute and write descriptors without fitting (`ml_type=-1`) |
| 2 | [`lml_fe_afs`](examples/lml_fe_afs) | LML fit, bcc Fe, AFS descriptor |
| 3 | [`qnml_hea_bso4`](examples/qnml_hea_bso4) | QNML fit, Ta-Ti-V-W HEA, bispectrum SO(4) |
| 4 | [`kernel_poly_w_bso4`](examples/kernel_poly_w_bso4) | Kernel (KNML) fit, W, polynomial kernel |
| 5 | [`kernel_r_fe_bso4`](examples/kernel_r_fe_bso4) | Kernel fit with a random kernel, done in a single step |
| 6 | [`lml_afs_read_write_dump`](examples/lml_afs_read_write_dump) | Dump descriptors, then read them back for a fit |
| 7 | [`lml_fe_bso4_zbl`](examples/lml_fe_bso4_zbl) | LML fit, bcc Fe, with k2b and a ZBL repulsive core |
| 8 | [`lml_hea_kace`](examples/lml_hea_kace) | LML fit, HEA, kACE descriptor |
| 9 | [`qnml_hea_ace`](examples/qnml_hea_ace) | QNML fit, HEA, ACE descriptor |

## Tests

[`tests/`](tests) holds 26 regression fits covering LML, QNML, polynomial
chaos, kernels and descriptor writing, with the AFS, SO(3)/SO(4), PiP,
Fourier N-body and ACE/kACE descriptors, on LAPACK and ScaLapack. The list is
in [`tests/TESTS_SUMMARY.md`](tests/TESTS_SUMMARY.md).

```bash
cd tests
# first set mld_exe (path to milady_main.exe) and nprocs in local_small_test.py
python local_small_test.py run     # run all fits and compare them with REF/
python local_small_test.py clean   # remove generated files
```

## Repository layout

```text
src/            Fortran sources (mld/, ndm/, scatools/, …) and build scripts (src/scripts/)
cmake_files/    CMake modules
examples/       Step-by-step example runs
tests/          Regression test suite
tools/          Potential-validation tools for LAMMPS
```text

## Publications

If you use MILADY, please cite:

- A. M. Goryaeva, J. Dérès, C. Lapointe, P. Grigorev, T. D. Swinburne,
  J. R. Kermode, L. Ventelon, J. Baima, M.-C. Marinica, *Efficient and
  transferable machine learning potentials for the simulation of crystal
  defects in bcc Fe and W*, **Phys. Rev. Mater. 5, 103803 (2021)**.
- A. M. Goryaeva, J.-B. Maillet, M.-C. Marinica, *Towards better efficiency of
  interatomic linear machine learning potentials*, **Comp. Mater. Sci. 166,
  200 (2019)**.


## Contributors

MILADY development started in 2015 at SRMP, CEA Saclay, France. The current
architecture was designed by M.-C. Marinica and A. M. Goryaeva.

**Main contributors**, all current or former members of SRMP, CEA Saclay:
M.-C. Marinica (2015–present), A. Zhong (2022–present), C. Lapointe
(2018–present), W. Unn-Toc (2015–2017), A. M. Goryaeva (2018–2024), J. Deres
(2019–2021), J. Baima (2020–2022), A. Allera (2022–2024).

**Software development support:** Anida Khizar and Christian Van Wambeke
(LGLS, CEA Saclay).
**MILADY–LAMMPS coupling:** Thomas D. Swinburne (CINaM, Marseille).

Contributions are welcome. Please open an issue or a pull request.

## Licence

MILADY is distributed under the **Academic Software Licence (ASL)**. See
[`LICENSE.md`](LICENSE.md). The ASL is a source-available licence modelled on
GPLv2 that allows **academic, non-commercial use only**. For commercial use
rights, contact the rights holder.
