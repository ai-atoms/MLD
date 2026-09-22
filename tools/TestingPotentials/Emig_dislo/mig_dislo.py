import os, re, argparse
import shutil
import subprocess
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

from ase import Atoms
from ase.io import read, write
from typing import Tuple
from scipy.interpolate import CubicSpline

class Migration : 
    def __init__(self, 
                 a0 : float,
                 a0_ref : float,
                 species : str,
                 partition_lammps : str,
                 script_lammps : str,
                 mpi_command : str,
                 hostfile_content : str, 
                 path_lmp : os.PathLike[str],
                 path_pot : os.PathLike[str],
                 extra_pot_files : list = None) -> None : 
        
        self.a0 = a0
        self.a0_ref = a0_ref
        self.species = species
        self.partition_lammps = partition_lammps
        self.hostfile_content = hostfile_content
        self.script_lammps_min = f"""
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
        
        self.script_lammps_min +="""
            print "Number_of_atoms = ${natoms}"
            print "Energy_box = ${tenergy}"
            print "sig_xx = ${sxx}"
            print "sig_yy = ${syy}"
            print "sig_zz = ${szz}"
            print "sig_xy = ${sxy}"
            print "sig_xz = ${sxz}"
            print "sig_yz = ${syz}"
            print "Volume = ${vol_box}"

            write_data  relaxed.lmp"""

        self.script_lammps_neb = f"""
            units           metal
            boundary        p p p
            atom_style      atomic
            atom_modify     map array sort 0 0.0
            read_data        initial.lmp
            {script_lammps}

            minimize                 0.0 1.0e-4 1000 10000
            reset_timestep           0
            timestep        0.01 

            fix             1 all neb 15.0
            min_style       fire
            thermo_style custom step etotal pxx pyy pzz pxy pxz pyz
            thermo 5
            neb             0.0 1.0e-2 100 100 10 final final.data

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

        self.script_lammps_neb +="""
            print "Number_of_atoms = ${natoms}"
            print "Energy_box = ${tenergy}"
            print "sig_xx = ${sxx}"
            print "sig_yy = ${syy}"
            print "sig_zz = ${szz}"
            print "sig_xy = ${sxy}"
            print "sig_xz = ${sxz}"
            print "sig_yz = ${syz}"
            print "Volume = ${vol_box}"

            variable i equal part
            write_data  final_neb.$i"""    

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
        # Auto-discover migration paths: each subdirectory in config/
        config_dir = f'{self.init_dir}/config'
        self.path_mig = sorted([
            f'{config_dir}/{d}' for d in os.listdir(config_dir)
            if os.path.isdir(f'{config_dir}/{d}')
        ])

    def copy_extra_files(self, path_work : os.PathLike[str]) -> None :
        for f in self.extra_pot_files:
            if not os.path.exists(f):
                print(f"WARNING: extra potential file '{f}' not found, skipping")
                continue
            dest = f'{path_work}/{os.path.basename(f)}'
            if not os.path.exists(dest):
                os.symlink(f, dest)

    def clean_directories(self) -> None : 
        # Derive work directories from discovered migration configs
        dirs_to_delete = [
            os.path.basename(p) for p in self.path_mig
        ]
        for directory in dirs_to_delete : 
            try :
                shutil.rmtree(f'{self.init_dir}/{directory}')
                os.remove(f'{self.init_dir}/{directory}.png')
            except : 
                continue

        if os.path.exists(f'{self.init_dir}/mig_dislo.log') :
            os.remove(f'{self.init_dir}/mig_dislo.log')
                
        return

    def get_energy_from_lammps(self, out_lammps : os.PathLike[str]) -> float :
        with open(out_lammps, 'r') as file:
            for line in file:
                if line.startswith('Energy_box'):
                    return float(line.split('=')[1])
        return None

    def logger(self, name_log : str, txt : str) -> None : 
        with open(name_log, 'a') as w :
            w.write(f'{txt} \n')
            w.close()

        print(txt)
        return

    def dump_inputs_lammps(self, path : os.PathLike[str], mode : str = 'min') -> None: 
        with open(f'{path}/in.lammps','w') as w : 
            if mode == 'min' :
                w.write(self.script_lammps_min)
            if mode == 'neb' : 
                w.write(self.script_lammps_neb)
            w.close()
        return 

    def dump_final_neb_file(self, path_relax : os.PathLike[str], path_final : os.PathLike[str]) -> None : 
        """Parse a LAMMPS data file (orthogonal or triclinic) and write a NEB final-state file."""
        r = open(path_relax,'r').readlines() 
        nb_atoms = None
        cell = np.zeros((3,3))
        in_atoms = False
        compt_pos = 0

        # First pass: extract nb_atoms and cell bounds
        for l in r:
            tokens = l.split()
            if len(tokens) >= 2 and tokens[1] == 'atoms':
                nb_atoms = int(tokens[0])
            elif len(tokens) >= 4 and tokens[2] == 'xlo' and tokens[3] == 'xhi':
                cell[0,0] = float(tokens[1]) - float(tokens[0])
            elif len(tokens) >= 4 and tokens[2] == 'ylo' and tokens[3] == 'yhi':
                cell[1,1] = float(tokens[1]) - float(tokens[0])
            elif len(tokens) >= 4 and tokens[2] == 'zlo' and tokens[3] == 'zhi':
                cell[2,2] = float(tokens[1]) - float(tokens[0])

        with open(path_final, 'w') as w :
            w.write(f'{nb_atoms} \n')
            for l in r :
                tokens = l.split()
                if len(tokens) >= 1 and tokens[0] == 'Atoms':
                    in_atoms = True
                    continue
                if in_atoms and len(tokens) == 0:
                    continue  # skip blank line after "Atoms" header
                if in_atoms and len(tokens) >= 5:
                    # atom-id type x y z [ix iy iz]
                    # Write wrapped coordinates only (no image flag unwrapping)
                    # LAMMPS NEB handles periodicity internally
                    w.write(f'{int(tokens[0])} {tokens[2]} {tokens[3]} {tokens[4]} \n')
                    compt_pos += 1
                    if compt_pos >= nb_atoms:
                        break
        return      

    def scale_atoms(self, atoms : Atoms) -> Atoms : 
        """Uniform rescaling from the reference lattice parameter to the target one.
        Works for arbitrary (non-cubic) cells like dislocation dipoles."""
        ratio = self.a0 / self.a0_ref
        new_cell = atoms.cell[:] * ratio
        atoms.set_cell(new_cell, scale_atoms=True)
        return atoms
    
    def OneMinimisation(self, path_poscar : os.PathLike[str], path_mig : os.PathLike[str]) -> Tuple[float,int] :
        self.logger(f'{self.init_dir}/mig_dislo.log',
                    f'... Computing for {path_poscar.split("/")[-2]}')
        name_poscar = os.path.basename(path_poscar)
        path_work = f"{path_mig}/{name_poscar.split('.')[0]}"
        if not os.path.exists(path_work) : 
            os.mkdir(path_work)
        link_path = f'{path_work}/{self.pot_local_name}'
        if not os.path.exists(link_path):
            os.symlink(self.path_pot, link_path)
        self.copy_extra_files(path_work)

        # copying and writing input lammps ! 
        atoms_poscar = read(path_poscar)
        atoms_poscar = self.scale_atoms(atoms_poscar)
        atoms_poscar.set_chemical_symbols(len(atoms_poscar)*[self.species])
        write(f'{path_work}/in.lmp', atoms_poscar, format='lammps-data')
        self.dump_inputs_lammps(path_work)
        
        if self.hostfile_content:
            with open(f'{path_work}/hostfile', "w") as file:
                file.write(self.hostfile_content)
        
        # lammps time ! 
        lammps_command = f'{self.mpi_command} {self.path_lmp} < in.lammps > out.lammps'
        _ = subprocess.call(lammps_command,
                        shell=True,
                        cwd=path_work)
        energy_lammps = self.get_energy_from_lammps(f'{path_work}/out.lammps')
        self.logger(f'{self.init_dir}/mig_dislo.log',
                    f'... Energy {path_poscar.split("/")[-2]} = {energy_lammps} eV ')
        return energy_lammps, len(atoms_poscar)
    

    def Minimisation(self, path_migration : os.PathLike[str]) -> os.PathLike[str] : 
        name_migration = os.path.basename(path_migration)        
        path_work_mig = f'{self.init_dir}/{name_migration}'
        if not os.path.exists(path_work_mig) : 
            os.mkdir(path_work_mig)

        for mig_poscar in sorted([f'{path_migration}/{f}' for f in os.listdir(path_migration)]) : 
            if self.hostfile_content:
                with open(f'{path_work_mig}/hostfile', "w") as file:
                    file.write(self.hostfile_content)
            _, _ = self.OneMinimisation(mig_poscar, path_work_mig)
            self.logger(f'{self.init_dir}/mig_dislo.log',
                        f'---------------------------------------------------------------------- \n')      
        return path_work_mig
    

    def PerformNEB(self, path_migration : os.PathLike[str]) -> None : 
        #link pot
        link_path = f'{path_migration}/{self.pot_local_name}'
        if not os.path.exists(link_path):
            os.symlink(self.path_pot, link_path)
        self.copy_extra_files(path_migration)
        
        # put initial config (ini.poscar -> ini/ work dir -> relaxed.lmp)
        os.system(f'cp {path_migration}/ini/relaxed.lmp {path_migration}/initial.lmp')
        
        # put final config (end.poscar -> end/ work dir -> relaxed.lmp)
        self.dump_final_neb_file(f'{path_migration}/end/relaxed.lmp',
                                 f'{path_migration}/final.data')
        self.dump_inputs_lammps(path_migration, mode='neb')
        
        if self.hostfile_content:
            with open(f'{path_migration}/hostfile', "w") as file:
                file.write(self.hostfile_content)

        lammps_command = f'{self.mpi_command} {self.path_lmp} -partition {self.partition_lammps} -i in.lammps > out_neb.lammps'
        _ = subprocess.call(lammps_command,
                        shell=True,
                        cwd=path_migration)       
        return 
        
    def get_climbing_replica(self, out_neb : os.PathLike[str]) -> int : 
        with open(out_neb,'r') as r : 
            for line in r:
                if re.match('(.*)Climbing replica(.*)', line):
                    word=line.split()
                    return int(word[-1])

        return None

    def PlotBarrier(self, path_migration : os.PathLike[str]) -> None :
        list_screen = [f for f in os.listdir(path_migration) if f.startswith('screen.')]
        list_screen.sort(key=lambda f: int(f.split('.')[1]))
        list_coord, list_energy = [], []
        for id_screen, screen in enumerate(list_screen) :
            list_coord.append(float(id_screen+1))
            list_energy.append(self.get_energy_from_lammps(f'{path_migration}/{screen}'))

        coord = np.array(list_coord)/np.amax(list_coord)
        energy = np.array(list_energy)

        ref_energy = np.amin(energy)
        energy += - ref_energy
        cmin=coord[0]-0.05/100.0
        cmax=coord[-1]+0.05/100.0
        xepsilon=np.linspace(cmin,cmax,num=200,endpoint=True)
        y_atomic_smooth=CubicSpline(coord,energy)
        fig=plt.figure()
        ax=fig.add_subplot(111)
        ax.plot(coord,energy,'o',markersize=6,color='slateblue',label='atomic') 
        ax.plot(xepsilon,y_atomic_smooth(xepsilon),linestyle='--',markersize=6,color='slateblue', label='interpolation') 
        ax.set_title('NEB energy reference {:2.3f} eV'.format(ref_energy))
        ax.set_ylabel(r" $E (r)$ (eV) ", color='slateblue',fontsize='20')
        ax.set_xlabel(r" r (reaction coordinate) ")
        ax.set_xlim([cmin,cmax])
        ax.tick_params(axis='x',labelsize='16')
        ax.tick_params(axis='y',labelsize='16')
        ax.legend(loc='best') 
        defsize=fig.get_size_inches()
        fig.set_size_inches( (defsize[0]*1.5,defsize[1]*1.5) )
        fig.savefig(f'{self.init_dir}/{os.path.basename(path_migration)}.png')
        return  

    def OneGetEnergyBarrier(self, path_migration : os.PathLike[str]) -> None : 
        idx_climb = self.get_climbing_replica(f'{path_migration}/out_neb.lammps')
        energy_ini = self.get_energy_from_lammps(f'{path_migration}/ini/out.lammps')
        energy_saddle = self.get_energy_from_lammps(f'{path_migration}/screen.{int(idx_climb)}')
        self.logger(f'{self.init_dir}/mig_dislo.log',
                    f'... Migration energy : {os.path.basename(path_migration)} = {energy_saddle - energy_ini} (saddle - ini) eV')
        self.logger(f'{self.init_dir}/mig_dislo.log',
                    f'---------------------------------------------------------------------- \n')
        return 

    def GetEnergyBarriers(self) -> None : 
        results = []
        for mig in self.path_mig : 
            path_work_mig = self.Minimisation(mig)
            self.PerformNEB(path_work_mig)
            self.PlotBarrier(path_work_mig)
            self.OneGetEnergyBarrier(path_work_mig)

            # Collect result for summary table
            name = os.path.basename(mig)
            idx_climb = self.get_climbing_replica(f'{path_work_mig}/out_neb.lammps')
            energy_ini = self.get_energy_from_lammps(f'{path_work_mig}/ini/out.lammps')
            energy_saddle = self.get_energy_from_lammps(f'{path_work_mig}/screen.{int(idx_climb)}')
            results.append((name, energy_saddle - energy_ini))

        # Print summary table
        header = f'{"Config":<25s} {"E_mig (eV)":>12s}'
        sep = '-' * len(header)
        self.logger(f'{self.init_dir}/mig_dislo.log', f'\n{sep}')
        self.logger(f'{self.init_dir}/mig_dislo.log', header)
        self.logger(f'{self.init_dir}/mig_dislo.log', sep)
        for name, e_mig in results:
            line = f'{name:<25s} {e_mig:>12.6f}'
            self.logger(f'{self.init_dir}/mig_dislo.log', line)
        self.logger(f'{self.init_dir}/mig_dislo.log', sep)

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

def read_mig_input(mig_input_path):
    mig_data = {}
    with open(mig_input_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                mig_data[key.strip()] = value.split('#')[0].strip()
    return mig_data

########################################
### MAIN
########################################
parser = argparse.ArgumentParser('MigDislo')
parser.add_argument('-m','--mode',default="migration")
args = parser.parse_args()
mode = args.mode

# Read input files
tester_data = read_tester_input("../tester.input")
mig_data    = read_mig_input("./mig_dislo.input")

a0               = float(mig_data["lattice_parameter"])
a0_ref           = float(mig_data.get("reference_lattice_parameter", "2.8553"))
species          = tester_data["species"]
mpi_command      = tester_data["neb_mpi_command"]
partition_lammps = tester_data["neb_partition_lammps"]
hostfile_content = tester_data.get("hostfile_content", "")
path_lmp         = tester_data["lammps_executable"]
path_pot         = tester_data["potential_file"]
lammps_script    = tester_data["lammps_script_block"]

print(f"   ------> lattice_parameter     = {a0}")
print(f"   ------> reference_lattice_par = {a0_ref}")
print(f"   ------> species               = {species}")
print(f"   ------> mpi_command           = {mpi_command}")
print(f"   ------> partition_lammps      = {partition_lammps}")
print(f"   ------> lammps_executable     = {path_lmp}")
print(f"   ------> potential_file        = {path_pot}")

# Extra potential files (optional)
extra_pot_files = []
for key in ['potential_file_ext1', 'potential_file_ext2']:
    val = tester_data.get(key, "").strip()
    if val:
        extra_pot_files.append(val)
if extra_pot_files:
    print(f"   ------> extra_pot_files       = {extra_pot_files}")

obj_migration = Migration(a0,
                          a0_ref,
                          species,
                          partition_lammps,
                          lammps_script,
                          mpi_command,
                          hostfile_content,
                          path_lmp,
                          path_pot,
                          extra_pot_files)
if mode == 'migration' :
    obj_migration.GetEnergyBarriers()

elif mode == 'clean' :
    obj_migration.clean_directories()

else :
    raise NotImplementedError(f'Mode : {mode} is not implemented')
