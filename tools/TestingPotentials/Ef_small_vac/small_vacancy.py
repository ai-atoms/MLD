import os, argparse
import shutil
import subprocess
import numpy as np

from ase import Atoms
from ase.io import read, write
from typing import Tuple

class Minimizer : 
    def __init__(self, 
                 a0 : float,
                 species : str,
                 script_lammps : str,
                 mpi_command : str,
                 path_lmp : os.PathLike[str],
                 path_pot : os.PathLike[str],
                 extra_pot_files : list = None) -> None : 
        
        self.a0 = a0
        self.species = species
        self.script_lammps = f"""
            units           metal
            boundary        p p p
            atom_style      atomic
            atom_modify     map array sort 0 0.0
            read_data        in.lmp
            {script_lammps}

            timestep        0.01
            minimize                 1.0e-6  1.0e-3 1000 10000
            thermo 1
            thermo_style custom step etotal  fmax pxx pyy pzz pxy pxz pyz pe

            variable vol_box equal vol
            variable natoms equal  atoms
            variable        sxx equal pxx
            variable        syy equal pyy
            variable        szz equal pzz
            variable        sxy equal pxy
            variable        sxz equal pxz
            variable        syz equal pyz
            variable tenergy equal  pe
            run 0 \n"""
        
        self.script_lammps +="""
            print "Number_of_atoms = ${natoms}"
            print "Energy_box = ${tenergy}"
            print "sig_xx = ${sxx}"
            print "sig_yy = ${syy}"
            print "sig_zz = ${szz}"
            print "sig_xy = ${sxy}"
            print "sig_xz = ${sxz}"
            print "sig_yz = ${syz}"
            print "Volume = ${vol_box}"

            write_data  relaxed.data"""
        
        self.mpi_command = mpi_command
        self.path_lmp = path_lmp
        self.path_pot = path_pot
        self.extra_pot_files = extra_pot_files or []
        # Extract the local potential filename from pair_coeff line in the LAMMPS script
        self.pot_local_name = os.path.basename(path_pot)
        for line in script_lammps.splitlines():
            if 'pair_coeff' in line:
                tokens = line.split()
                idx = len(tokens) - 1
                while idx >= 0 and not tokens[idx].startswith('*'):
                    idx -= 1
                if idx + 1 < len(tokens):
                    self.pot_local_name = tokens[idx + 1]
                break
    
        self.init_dir = os.getcwd()
        self.path_poscar_bulk = f'{self.init_dir}/config/bulk.poscar'
        # Auto-discover defect POSCARs: everything in config/ except bulk.poscar
        config_dir = f'{self.init_dir}/config'
        self.path_poscar_dfct = sorted([
            f'{config_dir}/{f}' for f in os.listdir(config_dir)
            if f.endswith('.poscar') and f != 'bulk.poscar'
        ])
        self.reference_energy = None
        self.reference_nb_atom = None

    def copy_extra_files(self, path_work : os.PathLike[str]) -> None :
        for f in self.extra_pot_files:
            if not os.path.exists(f):
                print(f"WARNING: extra potential file '{f}' not found, skipping")
                continue
            dest = f'{path_work}/{os.path.basename(f)}'
            if not os.path.exists(dest):
                os.symlink(f, dest)

    def clean_directories(self) -> None : 
        # Derive work directories from the POSCAR filenames
        dirs_to_delete = ['bulk'] + [
            os.path.basename(p).split('.')[0] for p in self.path_poscar_dfct
        ]
        for directory in dirs_to_delete : 
            try :
                shutil.rmtree(f'{self.init_dir}/{directory}')
            except : 
                continue
        if os.path.exists(f'{self.init_dir}/small_vac.log') :
            os.remove(f'{self.init_dir}/small_vac.log')
        return 

    def get_energy_from_lammps(self, out_lammps : os.PathLike[str]) -> float :
        with open(out_lammps, 'r') as file:
            for line in file:
                if line.startswith('Energy_box'):
                    return float(line.split('=')[1])

    def logger(self, name_log : str, txt : str) -> None : 
        with open(name_log, 'a') as w :
            w.write(f'{txt} \n')
            w.close()

        print(txt)
        return

    def dump_inputs_lammps(self, path : os.PathLike[str]) -> None: 
        with open(f'{path}/in.lammps','w') as w : 
            w.write(self.script_lammps)
            w.close()
        return 

    def get_cell_size(self, atoms : Atoms) -> int :     
        return np.round( np.power((len(atoms)/2.0), 1./3.) )
    
    def scale_atoms(self, atoms : Atoms, size : int) -> Atoms : 
        new_cell = size*self.a0*np.eye(3)
        atoms.set_cell(new_cell, scale_atoms=True)
        return atoms
    
    def OneMinimisation(self, path_poscar : os.PathLike[str]) -> Tuple[float,int] :
        self.logger(f'{self.init_dir}/small_vac.log',
                    f'... Computing for {path_poscar}')
        name_poscar = os.path.basename(path_poscar)
        path_work = f"{self.init_dir}/{name_poscar.split('.')[0]}"
        if not os.path.exists(path_work) : 
            os.mkdir(path_work)
        
        link_path = f'{path_work}/{self.pot_local_name}'
        if not os.path.exists(link_path):
            os.symlink(self.path_pot, link_path)
        self.copy_extra_files(path_work)
        # copying and writing input lammps ! 
        atoms_poscar = read(path_poscar)
        size_system = self.get_cell_size(atoms_poscar)
        atoms_poscar = self.scale_atoms(atoms_poscar, size_system)
        atoms_poscar.set_chemical_symbols(len(atoms_poscar)*[self.species])
        write(f'{path_work}/in.lmp', atoms_poscar, format='lammps-data')
        self.dump_inputs_lammps(path_work)
        
        # lammps time ! 
        lammps_command = f'{self.mpi_command} {self.path_lmp} < in.lammps > out.lammps'
        _ = subprocess.call(lammps_command,
                        shell=True,
                        cwd=path_work)
        energy_lammps = self.get_energy_from_lammps(f'{path_work}/out.lammps')
        self.logger(f'{self.init_dir}/small_vac.log',
                    f'... Energy {path_poscar} = {energy_lammps} eV ')
        return energy_lammps, len(atoms_poscar)
    
    def formation_energy(self, e_dfct : float, n_dfct : int) -> float :
        return e_dfct - (float(n_dfct)/float(self.reference_nb_atom))*self.reference_energy

    def GetFormationEnergy(self) -> None : 
        self.reference_energy, self.reference_nb_atom = self.OneMinimisation(self.path_poscar_bulk)
        self.logger(f'{self.init_dir}/small_vac.log',
                    f'---------------------------------------------------------------------- \n')      

        # Collect results: list of (name, n_vac, E_f, E_bind)
        results = []
        ef_1vac = None

        for dfct_poscar in self.path_poscar_dfct : 
            dfct_energy, nb_at_dfct = self.OneMinimisation(dfct_poscar)
            ef = self.formation_energy(dfct_energy, nb_at_dfct)
            n_vac = self.reference_nb_atom - nb_at_dfct
            name = os.path.basename(dfct_poscar).replace('.poscar', '')

            if ef_1vac is None:
                ef_1vac = ef  # first defect is assumed to be 1vac

            e_bind = n_vac * ef_1vac - ef if n_vac > 1 else 0.0
            results.append((name, n_vac, ef, e_bind))

            self.logger(f'{self.init_dir}/small_vac.log',
                        f'... Formation energy : {name} = {ef} eV')
            self.logger(f'{self.init_dir}/small_vac.log',
                        f'---------------------------------------------------------------------- \n')      

        # Print summary table
        header = f'{"Config":<20s} {"N_vac":>5s} {"E_f (eV)":>12s} {"E_bind (eV)":>12s}'
        sep = '-' * len(header)
        self.logger(f'{self.init_dir}/small_vac.log', f'\n{sep}')
        self.logger(f'{self.init_dir}/small_vac.log', header)
        self.logger(f'{self.init_dir}/small_vac.log', sep)
        for name, n_vac, ef, e_bind in results:
            line = f'{name:<20s} {n_vac:>5d} {ef:>12.6f} {e_bind:>12.6f}'
            self.logger(f'{self.init_dir}/small_vac.log', line)
        self.logger(f'{self.init_dir}/small_vac.log', sep)

        return


########################################
### INPUT READERS
########################################
def read_tester_input(tester_input_path):
    tester_data = {}
    current_key = None
    current_value = []
    with open(tester_input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '"""' in line:
                if current_key:
                    tester_data[current_key] = '\n'.join(current_value).strip()
                    current_key = None
                    current_value = []
                else:
                    key, _ = line.split('=', 1)
                    current_key = key.strip()
                continue
            elif '=' in line and current_key is None:
                key, value = line.split('=', 1)
                tester_data[key.strip()] = value.strip()
            elif current_key:
                current_value.append(line)
    return tester_data

def read_local_input(local_input_path):
    local_data = {}
    with open(local_input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                local_data[key.strip()] = value.split('#')[0].strip()
    return local_data

########################################
### MAIN
########################################
parser = argparse.ArgumentParser('SmallVac')
parser.add_argument('-m','--mode',default="formation")
args = parser.parse_args()
mode = args.mode

# Read input files
tester_data = read_tester_input("../tester.input")
local_data  = read_local_input("./small_vacancy.input")

a0               = float(local_data["lattice_parameter"])
species          = tester_data["species"]
mpi_command      = tester_data["min_mpi_command"]
path_lmp         = tester_data["lammps_executable"]
path_pot         = tester_data["potential_file"]
lammps_script    = tester_data["lammps_script_block"]

print(f"   ------> lattice_parameter = {a0}")
print(f"   ------> species           = {species}")
print(f"   ------> mpi_command       = {mpi_command}")
print(f"   ------> lammps_executable = {path_lmp}")
print(f"   ------> potential_file    = {path_pot}")

# Extra potential files (optional)
extra_pot_files = []
for key in ['potential_file_ext1', 'potential_file_ext2']:
    val = tester_data.get(key, "").strip()
    if val:
        extra_pot_files.append(val)
if extra_pot_files:
    print(f"   ------> extra_pot_files   = {extra_pot_files}")

obj_minimizer = Minimizer(a0,
                          species,
                          lammps_script,
                          mpi_command,
                          path_lmp,
                          path_pot,
                          extra_pot_files)
if mode == 'formation' :
    obj_minimizer.GetFormationEnergy()

elif mode == 'clean' :
    obj_minimizer.clean_directories()

else :
    raise NotImplementedError(f'Mode : {mode} is not implemented')
