#!/usr/bin/env python
"""
MiLaDy fit test runner with filtering capabilities.

Usage:
  python run_tests.py <mode> [--exe <path>] [--nprocs <N>] [--model MODEL] [--descriptor DESC]
                      [--diag DIAG] [--element ELEM] [--subtype SUB]

Required arguments:
  mode              'run' to execute and check tests
                    'build' to run and save outputs as new references
                    'clean' to clean directories
  --exe   EXE       path to milady_main.exe  (default: ../mld_build_mix/bin/milady_main.exe)
  --nprocs N        number of MPI processes   (default: 6)

Optional filters (if not given, all values are included):
  --model       filter by model type:  LML, KER, QNML, PolyC, WriteDesc
  --descriptor  filter by descriptor:  BSO4, AFS, PSO3, PSO4, PiP, ACE, FTNBODY, ...
  --diag        filter by diag method: Lapack, Scalapack
  --element     filter by element:     W, HEA
  --subtype     filter by subtype:     tradACE, kACE, 3chBSO4, ...

Examples:
  python run_tests.py run
  python run_tests.py run --exe /path/to/milady_main.exe --nprocs 4
  python run_tests.py run --model LML
  python run_tests.py run --descriptor ACE --element W
  python run_tests.py run --model KER --diag Scalapack
  python run_tests.py build --descriptor ACE
  python run_tests.py clean
"""

import os
import re
import sys
import argparse
import glob
import shutil
import numpy as np
import subprocess


# ----------------------------------------------------------------
# Test registry: each entry holds the metadata and run/check logic
# ----------------------------------------------------------------
TEST_REGISTRY = {
    'fit01':  {'model': 'LML',       'element': 'W',   'descriptor': 'BS04',    'subtype': 'hybridG2',       'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit02':  {'model': 'LML',       'element': 'W',   'descriptor': 'AFS',     'subtype': 'AFS',            'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit03':  {'model': 'QNML',      'element': 'W',   'descriptor': 'BSO4',    'subtype': 'BSO4',           'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit04':  {'model': 'QNML',      'element': 'W',   'descriptor': 'BSO4',    'subtype': 'BSO4',           'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit05':  {'model': 'PolyC',     'element': 'W',   'descriptor': 'BSO4',    'subtype': 'BSO4',           'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit06':  {'model': 'PolyC',     'element': 'W',   'descriptor': 'BSO4',    'subtype': 'BSO4',           'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit07':  {'model': 'QNML',      'element': 'HEA', 'descriptor': 'PSO4',    'subtype': '3chPSO4',        'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit08':  {'model': 'LML',       'element': 'HEA', 'descriptor': 'BSO4',    'subtype': '3chBSO4',        'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit09':  {'model': 'LML',       'element': 'HEA', 'descriptor': 'AFS',     'subtype': '3chAFS',         'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit10':  {'model': 'KER',       'element': 'HEA', 'descriptor': 'PSO4',    'subtype': '3chPSO4',        'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit11':  {'model': 'LML',       'element': 'W',   'descriptor': 'PSO3',    'subtype': 'PSO3',           'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit12':  {'model': 'LML',       'element': 'W',   'descriptor': 'AFS',     'subtype': 'WRITE/READ DUMP','diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': 'dump'},
    'fit13':  {'model': 'KER',       'element': 'W',   'descriptor': 'BSO4',    'subtype': 'BSO4',           'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': '3step'},
    'fit14':  {'model': 'KER',       'element': 'W',   'descriptor': 'BSO4',    'subtype': 'BSO4',           'diag': 'Scalapack', 'check': 'fit',       'desc_cols': None, 'multistep': '3step'},
    'fit15':  {'model': 'KER',       'element': 'W',   'descriptor': 'BSO4',    'subtype': 'BSO4',           'diag': 'Scalapack', 'check': 'fit',       'desc_cols': None, 'multistep': '3step'},
    'fit16':  {'model': 'LML',       'element': 'W',   'descriptor': 'PiP',     'subtype': 'PiP',            'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit17':  {'model': 'WriteDesc', 'element': 'W',   'descriptor': 'BSO4',    'subtype': 'BSO4',           'diag': 'Lapack',    'check': 'desc',      'desc_cols': 56,   'multistep': False},
    'fit18':  {'model': 'WriteDesc', 'element': 'W',   'descriptor': 'AFS',     'subtype': 'AFS',            'diag': 'Lapack',    'check': 'desc',      'desc_cols': 133,  'multistep': False},
    'fit19':  {'model': 'WriteDesc', 'element': 'HEA', 'descriptor': 'BSO4',    'subtype': '3chBSO4',        'diag': 'Lapack',    'check': 'desc',      'desc_cols': 166,  'multistep': False},
    'fit20':  {'model': 'LML',       'element': 'W',   'descriptor': 'FTNBODY', 'subtype': 'FTNBODY',        'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit21':  {'model': 'LML',       'element': 'W',   'descriptor': 'ACE',     'subtype': 'tradACE',        'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit22':  {'model': 'LML',       'element': 'W',   'descriptor': 'ACE',     'subtype': 'kACE',           'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit23':  {'model': 'LML',       'element': 'W',   'descriptor': 'ACE',     'subtype': 'k2btradACE',     'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit24':  {'model': 'LML',       'element': 'HEA', 'descriptor': 'ACE',     'subtype': 'tradACE',        'diag': 'Lapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
    'fit25':  {'model': 'LML',       'element': 'HEA', 'descriptor': 'ACE',     'subtype': 'kACE',           'diag': 'Scalapack',    'check': 'fit',       'desc_cols': None, 'multistep': False},
}

# ----------------------------------------------------------------
# Helper functions (same logic as local_small_test.py)
# ----------------------------------------------------------------

def run_ndm(mld_exe, nprocs):
    command = ["mpirun", "-np", str(nprocs), str(mld_exe)]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode == 0:
        with open('ndm.out', 'w') as f:
            f.write(result.stdout.decode())
    else:
        print("    Execution failed.")
        print("    Error:", result.stderr.decode())


def read_x(filename, i):
    return np.loadtxt(filename, usecols=(i))


def check_path(path):
    if os.path.isdir(path) and not os.path.islink(path):
        return 1
    elif os.path.islink(path):
        return 2
    else:
        return 0


def process_files(directory_path, desired_columns):
    statuss = check_path(directory_path)
    if statuss == 0:
        print(f"    {directory_path} does not exist")
        return np.array([])
    all_files = os.listdir(directory_path)
    eml_files = [f for f in all_files if f.endswith('.eml')]

    content_dict = {}
    for fname in eml_files:
        file_path = os.path.join(directory_path, fname)
        with open(file_path, 'r') as f:
            content = [line.split() for line in f]
            if len(content[0]) == desired_columns:
                num_lines = len(content)
                if num_lines not in content_dict:
                    content_dict[num_lines] = []
                content_dict[num_lines].append(content)
            else:
                print(f"    File {fname} has {len(content[0])} columns instead of {desired_columns}")

    results = {}
    for num_lines, matrices in content_dict.items():
        concatenated = np.vstack(matrices)
        norms = np.linalg.norm(concatenated[:, 1:], axis=1)
        results[num_lines] = np.mean(norms)

    sorted_results = sorted(results.items(), key=lambda x: x[0])
    return np.array(sorted_results)


class Fit(object):
    def __init__(self):
        self._energy_train = []
        self._force_train = []
        self._stress_train = []
        self._energy_dft = []
        self._force_dft = []
        self._stress_dft = []

    @property
    def energy_dft(self):
        return self._energy_dft
    @energy_dft.setter
    def energy_dft(self, fileE):
        self._energy_dft = read_x(fileE, 1)

    @property
    def energy_train(self):
        return self._energy_train
    @energy_train.setter
    def energy_train(self, fileE):
        self._energy_train = read_x(fileE, 0)

    @property
    def force_dft(self):
        return self._force_dft
    @force_dft.setter
    def force_dft(self, fileE):
        self._force_dft = read_x(fileE, 1)

    @property
    def force_train(self):
        return self._force_train
    @force_train.setter
    def force_train(self, fileE):
        self._force_train = read_x(fileE, 0)

    @property
    def stress_dft(self):
        return self._stress_dft
    @stress_dft.setter
    def stress_dft(self, fileE):
        self._stress_dft = read_x(fileE, 1)

    @property
    def stress_train(self):
        return self._stress_train
    @stress_train.setter
    def stress_train(self, fileE):
        self._stress_train = read_x(fileE, 0)


def compare_fit(fit_limit):
    cur_fit = Fit()
    ref_fit = Fit()

    cur_fit.energy_train = 'train_energy.out'
    ref_fit.energy_train = 'REF/train_energy.out'
    corr = np.corrcoef(cur_fit.energy_train, ref_fit.energy_train)[1, 0]
    print(f"    energy  correlation {corr:1.5f}", end="")
    if abs(corr - 1.0) < fit_limit:
        print("  OK")
    else:
        print("  FAILED")

    cur_fit.force_train = 'train_force.out'
    ref_fit.force_train = 'REF/train_force.out'
    corr = np.corrcoef(cur_fit.force_train, ref_fit.force_train)[1, 0]
    print(f"    force   correlation {corr:1.5f}", end="")
    if abs(corr - 1.0) < fit_limit:
        print("  OK")
    else:
        print("  FAILED")

    cur_fit.stress_train = 'train_stress.out'
    ref_fit.stress_train = 'REF/train_stress.out'
    corr = np.corrcoef(cur_fit.stress_train, ref_fit.stress_train)[1, 0]
    print(f"    stress  correlation {corr:1.5f}", end="")
    if abs(corr - 1.0) < fit_limit:
        print("  OK")
    else:
        print("  FAILED")


def compare_write_desc(vout, vref):
    for i in range(len(vout)):
        ff = 1e-10
        vdiff = abs(vout[i][1] - vref[i][1])
        if abs(vdiff) < ff:
            print(f"    diff {int(vout[i][0]):>5}   OK")
        else:
            print(f"    diff {int(vout[i][0]):>5}   FAILED  {vdiff:.18f}")


# ----------------------------------------------------------------
# Run / check logic for a single test
# ----------------------------------------------------------------

def execute_test(test_name, info, mld_exe, nprocs):
    """Run milady and check results for one test."""
    multistep = info['multistep']

    if multistep == 'dump':
        os.system("cp vacancy.ml_write_dump  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_read_dump  vacancy.ml")
        run_ndm(mld_exe, nprocs)
    elif multistep == '3step':
        os.system("cp vacancy.ml_01  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_02  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_03  vacancy.ml")
        run_ndm(mld_exe, nprocs)
    else:
        run_ndm(mld_exe, nprocs)

    # Validation
    if info['check'] == 'desc':
        vout = process_files('descDB/', info['desc_cols'])
        vref = process_files('REF/descDB/', info['desc_cols'])
        compare_write_desc(vout, vref)
    else:
        compare_fit(1.e-3)


def build_references(test_name, info, mld_exe, nprocs):
    """Run milady and copy outputs into REF/ as new reference data."""
    multistep = info['multistep']

    # Run milady (same logic as execute_test)
    if multistep == 'dump':
        os.system("cp vacancy.ml_write_dump  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_read_dump  vacancy.ml")
        run_ndm(mld_exe, nprocs)
    elif multistep == '3step':
        os.system("cp vacancy.ml_01  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_02  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_03  vacancy.ml")
        run_ndm(mld_exe, nprocs)
    else:
        run_ndm(mld_exe, nprocs)

    # Create REF/ if it does not exist
    os.makedirs('REF', exist_ok=True)

    # Copy *.out and *.milady files into REF/
    for pattern in ['*.out', '*.milady']:
        for f in glob.glob(pattern):
            shutil.copy2(f, 'REF/')
            print(f"    copied {f} -> REF/{f}")

    # Copy ndm.out
    if os.path.exists('ndm.out'):
        shutil.copy2('ndm.out', 'REF/')
        print("    copied ndm.out -> REF/ndm.out")

    # For WriteDesc tests, copy descDB/*.eml into REF/descDB/
    if info['check'] == 'desc':
        os.makedirs('REF/descDB', exist_ok=True)
        # Remove old reference eml files
        for old in glob.glob('REF/descDB/*.eml'):
            os.remove(old)
        for f in glob.glob('descDB/*.eml'):
            shutil.copy2(f, 'REF/descDB/')
        neml = len(glob.glob('REF/descDB/*.eml'))
        print(f"    copied {neml} .eml files -> REF/descDB/")

    print("    References built.")


# ----------------------------------------------------------------
# Filtering
# ----------------------------------------------------------------

def filter_tests(registry, model=None, descriptor=None, diag=None, element=None, subtype=None):
    """Return the subset of test names matching all given filters (case-insensitive)."""
    selected = []
    for name in sorted(registry.keys()):
        info = registry[name]
        if model and info['model'].lower() != model.lower():
            continue
        if descriptor and info['descriptor'].lower() != descriptor.lower():
            continue
        if diag and info['diag'].lower() != diag.lower():
            continue
        if element and info['element'].lower() != element.lower():
            continue
        if subtype and info['subtype'].lower() != subtype.lower():
            continue
        selected.append(name)
    return selected


# ----------------------------------------------------------------
# Main
# ----------------------------------------------------------------

def main():
    # Collect all valid values from registry for help messages
    all_models      = sorted(set(v['model']      for v in TEST_REGISTRY.values()))
    all_descriptors = sorted(set(v['descriptor']  for v in TEST_REGISTRY.values()))
    all_diags       = sorted(set(v['diag']        for v in TEST_REGISTRY.values()))
    all_elements    = sorted(set(v['element']     for v in TEST_REGISTRY.values()))
    all_subtypes    = sorted(set(v['subtype']     for v in TEST_REGISTRY.values()))

    parser = argparse.ArgumentParser(
        description="MiLaDy fit test runner with filtering.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('mode', choices=['run', 'build', 'clean'],
                        help="'run' to execute and check tests, 'build' to generate references, 'clean' to clean")
    parser.add_argument('--exe', default=None,
                        help="path to milady_main.exe (default: auto-detect from MLD_INSDIR or mld_build_mix)")
    parser.add_argument('--nprocs', type=int, default=6,
                        help="number of MPI processes (default: 6)")
    parser.add_argument('--model', default=None,
                        help=f"filter by model (all if omitted). Choices: {', '.join(all_models)}")
    parser.add_argument('--descriptor', default=None,
                        help=f"filter by descriptor (all if omitted). Choices: {', '.join(all_descriptors)}")
    parser.add_argument('--diag', default=None,
                        help=f"filter by diag (all if omitted). Choices: {', '.join(all_diags)}")
    parser.add_argument('--element', default=None,
                        help=f"filter by element (all if omitted). Choices: {', '.join(all_elements)}")
    parser.add_argument('--subtype', default=None,
                        help=f"filter by subtype (all if omitted). Choices: {', '.join(all_subtypes)}")

    args = parser.parse_args()

    # Resolve executable path
    if args.exe is not None:
        mld_exe = os.path.abspath(args.exe)
    elif 'MLD_INSDIR' in os.environ:
        mld_exe = os.path.join(os.environ['MLD_INSDIR'], 'bin', 'milady_main.exe')
    else:
        # Default relative to typical build location
        mld_exe = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'mld_build_mix', 'bin', 'milady_main.exe'))

    nprocs = args.nprocs
    BaseDir = os.getcwd()

    print(f"Executable: {mld_exe}")
    print(f"MPI procs:  {nprocs}")
    print()

    # Apply filters
    selected = filter_tests(TEST_REGISTRY,
                            model=args.model,
                            descriptor=args.descriptor,
                            diag=args.diag,
                            element=args.element,
                            subtype=args.subtype)

    if not selected:
        print("No tests match the given filters.")
        sys.exit(0)

    print(f"Selected {len(selected)} test(s): {', '.join(selected)}")
    print()

    # Extract databases before running
    if args.mode in ('run', 'build'):
        os.chdir(os.path.join(BaseDir, 'train_tests', 'data_for_tests'))
        os.system('tar -zxvf DB_HEA.tgz > /dev/null')
        os.system('tar -zxvf DB_one.tgz > /dev/null')
        os.system('tar -zxvf DB_two.tgz > /dev/null')

    if args.mode == 'clean':
        os.chdir(os.path.join(BaseDir, 'train_tests', 'data_for_tests'))
        os.system('rm -rf DB_HEA')
        os.system('rm -rf DB_one')
        os.system('rm -rf DB_two')

    # Run each selected test
    for i, test_name in enumerate(selected):
        info = TEST_REGISTRY[test_name]
        label = f"{info['model']:9s} {info['element']:3s} {info['descriptor']:7s} {info['subtype']}"
        print(f"[{i+1:2d}/{len(selected):2d}] {test_name}  {label}")

        test_dir = os.path.join(BaseDir, 'train_tests', test_name)
        os.chdir(test_dir)

        if args.mode == 'run':
            # Print README
            if os.path.exists('README'):
                with open('README', 'r') as f:
                    for line in f:
                        print(f"    {line.rstrip()}")

            execute_test(test_name, info, mld_exe, nprocs)

        if args.mode == 'build':
            # Print README
            if os.path.exists('README'):
                with open('README', 'r') as f:
                    for line in f:
                        print(f"    {line.rstrip()}")

            build_references(test_name, info, mld_exe, nprocs)

        if args.mode == 'clean':
            os.system('./clean.sh')

        print()

    print("Done.")


if __name__ == '__main__':
    main()
