# Summary of Fit Tests

| Test   | Model     | Element | Descriptor | Subtype         | Diag       |
|--------|-----------|---------|------------|-----------------|------------|
| fit01  | LML       | W       | BS04       | hybridG2        | Lapack     |
| fit02  | LML       | W       | AFS        | AFS             | Lapack     |
| fit03  | QNML      | W       | BSO4       | BSO4            | Lapack     |
| fit04  | QNML      | W       | BSO4       | BSO4            | Lapack     |
| fit05  | PolyC     | W       | BSO4       | BSO4            | Lapack     |
| fit06  | PolyC     | W       | BSO4       | BSO4            | Lapack     |
| fit07  | QNML      | HEA     | PSO4       | 3chPSO4         | Lapack     |
| fit08  | LML       | HEA     | BSO4       | 3chBSO4         | Lapack     |
| fit09  | LML       | HEA     | AFS        | 3chAFS          | Lapack     |
| fit10  | KER       | HEA     | PSO4       | 3chPSO4         | Lapack     |
| fit11  | LML       | W       | PSO3       | PSO3            | Lapack     |
| fit12  | LML       | W       | AFS        | WRITE/READ DUMP | Lapack     |
| fit13  | KER       | W       | BSO4       | BSO4            | Lapack     |
| fit14  | KER       | W       | BSO4       | BSO4            | Scalapack  |
| fit15  | KER       | W       | BSO4       | BSO4            | Scalapack  |
| fit16  | LML       | W       | PiP        | PiP             | Lapack     |
| fit17  | WriteDesc | W       | BSO4       | BSO4            | Lapack     |
| fit18  | WriteDesc | W       | AFS        | AFS             | Lapack     |
| fit19  | WriteDesc | HEA     | BSO4       | 3chBSO4         | Lapack     |
| fit20  | LML       | W       | FTNBODY    | FTNBODY         | Lapack     |
| fit21  | LML       | W       | ACE        | tradACE         | Lapack     |
| fit22  | LML       | W       | ACE        | kACE            | Lapack     |
| fit23  | LML       | W       | ACE        | k2btradACE      | Lapack     |
| fit24  | LML       | HEA     | ACE        | tradACE         | Lapack     |
| fit25  | LML       | HEA     | ACE        | kACE            | Lapack     |
| fit26  | LML       | HEA     | ACE        | k2bkACE         | ScaLapack  |

## How to Run the Tests

### Prerequisites

1. **Milady executable**: the path to `milady_main.exe` must be set in `local_small_test.py`
   (variable `mld_exe`). By default it points to:
   ```
   /home/marinica/MLD_GitHub/mld_build_mix/bin/milady_main.exe
   ```
   Adjust this path to your own build directory.

2. **MPI**: the tests are launched with `mpirun`. The number of MPI processes is set by
   `nprocs` in `local_small_test.py` (default: 6).

3. **Python**: requires `numpy` and standard library modules.

### Running the tests

From the `tests/` directory:

```bash
cd tests/
python local_small_test.py run
```

This will:
1. Extract the test databases (`DB_HEA.tgz`, `DB_one.tgz`, `DB_two.tgz`) in `train_tests/data_for_tests/`.
2. Enter each `train_tests/fitXX/` directory and run `milady_main.exe` with `mpirun -np 6`.
3. Compare outputs against reference data stored in each `fitXX/REF/` subdirectory.

### Cleaning

```bash
python local_small_test.py clean
```

This removes all generated output files and the extracted databases.

### Test structure

Each `fitXX/` directory contains:
- `vacancy.ml` — the MiLaDy input file (ML parameters).
- `vacancy.gin` — the geometry/NDM input file.
- `db_model.in` — the database configuration file.
- `name.in` — the name file.
- `README` — one-line description of the test.
- `REF/` — reference outputs (`train_energy.out`, `train_force.out`, `train_stress.out`, or `descDB/*.eml`).
- `clean.sh` — script to clean generated files.

### Validation criteria

- **Fit tests** (fit01–fit16, fit20–fit25): the correlation between current and reference
  `train_energy.out`, `train_force.out`, and `train_stress.out` is computed.
  The test passes if `|correlation - 1.0| < 1e-3`.
- **WriteDesc tests** (fit17–fit19): descriptor files in `descDB/` are compared against
  `REF/descDB/`. The norm of descriptor vectors is compared with a tolerance of `1e-10`.
- **Multi-step tests** (fit12–fit15): several `vacancy.ml_XX` files are run sequentially
  (e.g., write dump then read dump, or train/test/predict cycles).

---

## Advanced Test Runner: `run_tests.py`

The script `run_tests.py` is an improved version of `local_small_test.py` that supports
command-line arguments and filtering. It allows running a subset of tests based on model
type, descriptor, element, diagonalization method, or subtype.

### Usage

```bash
python run_tests.py <mode> [--exe EXE] [--nprocs N] [--model M] [--descriptor D]
                           [--diag G] [--element E] [--subtype S]
```

### Arguments

| Argument        | Required | Default                                | Description                                    |
|-----------------|----------|----------------------------------------|------------------------------------------------|
| `mode`          | yes      | —                                      | `run`, `build`, or `clean`                     |
| `--exe`         | no       | auto-detected from `MLD_INSDIR` or `../../mld_build_mix/bin/milady_main.exe` | Path to `milady_main.exe` |
| `--nprocs`      | no       | 6                                      | Number of MPI processes                        |
| `--model`       | no       | all                                    | Filter: LML, KER, QNML, PolyC, WriteDesc      |
| `--descriptor`  | no       | all                                    | Filter: BSO4, AFS, PSO3, PSO4, PiP, ACE, ...  |
| `--diag`        | no       | all                                    | Filter: Lapack, Scalapack                      |
| `--element`     | no       | all                                    | Filter: W, HEA                                |
| `--subtype`     | no       | all                                    | Filter: tradACE, kACE, 3chBSO4, ...           |

When a filter is **not specified**, all values for that category are included.
Multiple filters can be combined — only tests matching **all** given filters will run.

### Examples

```bash
# Run all tests with default executable and 6 MPI procs
python run_tests.py run

# Run with a custom executable and 4 procs
python run_tests.py run --exe /path/to/milady_main.exe --nprocs 4

# Run only LML model tests
python run_tests.py run --model LML

# Run only ACE descriptor tests on W
python run_tests.py run --descriptor ACE --element W

# Run only Scalapack diagonalization tests
python run_tests.py run --diag Scalapack

# Run only KER model with BSO4 descriptor
python run_tests.py run --model KER --descriptor BSO4

# Build new references for all ACE tests
python run_tests.py build --descriptor ACE

# Build new references for all tests (e.g. after changing -O flag)
python run_tests.py build

# Clean all test directories
python run_tests.py clean
```

### How it works

1. A **test registry** (Python dictionary) at the top of the script defines every test
   with its metadata: model, element, descriptor, subtype, diag method, validation type,
   and multi-step behaviour.
2. Command-line filters are applied to select the subset of tests to run.
3. Test databases are extracted from `train_tests/data_for_tests/*.tgz`.
4. For each selected test, the script:
   - Enters the `train_tests/fitXX/` directory.
   - Runs `mpirun -np <nprocs> <exe>` (possibly multiple times for multi-step tests).
   - **`run` mode**: compares outputs against reference data in `REF/`.
   - **`build` mode**: copies all outputs (`*.out`, `*.milady`, `ndm.out`, and `descDB/*.eml`
     for WriteDesc tests) into `REF/` as new references.
   - **`clean` mode**: runs `clean.sh` in each test directory.
5. Results are printed as OK or FAILED for each validation criterion (`run`),
   or as a list of copied files (`build`).

### Adding a new test

1. Create a new `train_tests/fitXX/` directory with the standard files.
2. Add a one-line entry to the `TEST_REGISTRY` dictionary in `run_tests.py`.
3. Update the summary table at the top of this file.


