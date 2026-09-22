import os, argparse 
import sys

# Function to read the tester.input file
def read_tester_input(tester_input_path):
    tester_data = {}
    current_key = None
    current_value = []

    with open(tester_input_path, 'r') as f:
        for line in f:
            line = line.strip()

            # Ignore comments and empty lines
            if not line or line.startswith('#'):
                continue

            # Detect the start of a multi-line block (starts with '"""' or '=')
            if '"""' in line:
                if current_key:  # We're already in a multi-line block, so close it
                    tester_data[current_key] = '\n'.join(current_value).strip()
                    current_key = None
                    current_value = []
                else:  # We're starting a new multi-line block
                    key, _ = line.split('=', 1)
                    current_key = key.strip()
                continue

            # Handle regular single-line key-value pairs
            elif '=' in line and current_key is None:
                key, value = line.split('=', 1)
                tester_data[key.strip()] = value.strip()

            # Handle lines that are part of a multi-line block
            elif current_key:
                current_value.append(line)

    #debug # Print what has been read for debugging
    #debug print("Read from tester.input:")
    #debug for key, value in tester_data.items():
    #debug     if '\n' in value:  # Detect multi-line values
    #debug         print(f"  {key}: \n{value}")
    #debug     else:
    #debug         print(f"  {key}: {value}")

    return tester_data



# Function to read the get_a0.input file
def read_get_a0_input(get_a0_input_path):
    get_a0_data = {}
    with open(get_a0_input_path, 'r') as f:
        for line in f:
            if '=' in line:
                key, value = line.split('=')
                get_a0_data[key.strip()] = value.strip()
    #debug print("\nRead from get_a0.input:")
    #debug for key, value in get_a0_data.items():
    #debug     print(f"  {key}: {value}")
    return get_a0_data


# Function to validate paths and inputs
def validate_inputs(tester_data, get_a0_data):
    # Check if lammps executable exists
    lammps_executable = tester_data.get("lammps_executable")
    if not lammps_executable or not os.path.exists(lammps_executable):
        print(f"Error: LAMMPS executable not found at {lammps_executable}")
        return False

    # Check if potential file exists
    potential_file = tester_data.get("potential_file")
    if not potential_file or not os.path.exists(potential_file):
        print(f"Error: Potential file not found at {potential_file}")
        return False

    # Check lattice parameter and remove comment if any
    lattice_param_raw = get_a0_data.get("lattice_parameter")
    #print(f"\n     -----> guess lattice parameter value = '{lattice_param_raw}'")  # Print raw value for debugging

    try:
        lattice_param = float(lattice_param_raw.split('#')[0].strip())  # Remove comment and convert to float
        if lattice_param <= 0:
            raise ValueError
    except (ValueError, TypeError):
        print(f"Error: Invalid lattice parameter. Received '{lattice_param_raw}'")
        return False

    print(f"   ------>   Lattice parameter read: {lattice_param}")

    # Check structure type
    structure_type = get_a0_data.get("structure_type").split('#')[0].strip()  # Remove comment
    print(f"   ------>   Structure type read: {structure_type}")
    valid_structures = ["bcc", "fcc", "a15", "c15"]
    
    # Split on "and" and validate each structure
    structure_list = [s.strip() for s in structure_type.split("and")]
    
    for struct in structure_list:
        if struct not in valid_structures:
            print(f"Error: Invalid structure type '{struct}'. Must be one of {valid_structures} or combinations using 'and'.")
            return False

    return True


# Function to write the LAMMPS input file
def create_lammps_input_file(tester_data, output_path):
    # Fixed parts of the LAMMPS script
    fixed_initial_part = """
clear
dimension 3
units		metal
boundary	p p p
atom_style	atomic
atom_modify     map array  sort 0 0.0 
box             tilt large
read_data	cube.lmp
neighbor 0.7 bin
neigh_modify delay 10 check yes 
"""

    fixed_final_part = """
min_style cg
fix 1 all nve
thermo_style custom step etotal pxx pyy pzz pxy pxz pyz
thermo 5
minimize                 0.0 1.0e-6 1000 10000
thermo			500
timestep	0.01
min_style	fire
variable vol_box equal vol
variable natoms equal  atoms
variable tenergy equal  pe
variable sxx equal pxx
variable syy equal pyy
variable szz equal pzz
variable sxy equal pxy
variable sxz equal pxz
variable syz equal pyz
print "Number_of_atoms = ${natoms}"
print "Energy_box = ${tenergy}"
print "sig_xx = ${sxx}"
print "sig_yy = ${syy}"
print "sig_zz = ${szz}"
print "sig_xy = ${sxy}"
print "sig_xz = ${sxz}"
print "sig_yz = ${syz}"
print "Volume = ${vol_box}"
write_data lammps.data
"""

    # Get the block from tester.input
    lammps_script_block = tester_data.get("lammps_script_block", "")
    
    # Ensure the directory exists
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Write the LAMMPS input file
    with open(output_path, 'w') as f:
        f.write(fixed_initial_part)
        f.write("\n")
        f.write(lammps_script_block)  # This now handles the multi-line block correctly
        f.write("\n")
        f.write(fixed_final_part)
    
    print(f"   ------>   LAMMPS input file created at {output_path}")


def modify_elastic_input(get_a0_data, tester_data, input_org_path, output_path):
    """
    Modify elastic/input_elastic_ini_org by replacing the placeholders with values
    from get_a0.input and tester.input and write the result to elastic/input_elastic_ini.

    - get_a0_data: Dictionary containing values from get_a0.input.
    - tester_data: Dictionary containing values from tester.input.
    - input_org_path: Path to the original input file (elastic/input_elastic_ini_org).
    - output_path: Path to write the modified input file (elastic/input_elastic_ini).
    """
    # Read the values from get_a0.input
    structure_type = get_a0_data.get("structure_type").split('#')[0].strip()  # Remove comment
    lattice_parameter = get_a0_data.get("lattice_parameter").split('#')[0].strip()  # Remove comment

    # Read the values from tester.input
    lammps_executable = tester_data.get("lammps_executable")
    potential_file = tester_data.get("potential_file")

    # Open the original file for reading and create a new file for writing the modified content
    with open(input_org_path, 'r') as file_org:
        content = file_org.read()

    # Replace the placeholders in the content
    content = content.replace("structure=MY_STRUCTURE", f"structure={structure_type}")
    content = content.replace("a_guess = MY_GUESS", f"a_guess = {lattice_parameter}")
    content = content.replace("lammps_exe= MY_EXE_LAMMPS", f"lammps_exe= {lammps_executable}")
    content = content.replace("input_potcar= MY_MLD_POT", f"input_potcar= {potential_file}")

    # Write the modified content to the output file
    with open(output_path, 'w') as file_output:
        file_output.write(content)

    print(f"   ------>   Modified elastic input file created at {output_path}")

import subprocess

def run_elastic_build():
    """
    Function to run the command `python Elastic.git/Src/elastic.py build`
    using subprocess.run.
    """
    command = ['python', 'Elastic.git/Src/elastic.py', 'build']
    os.chdir('./elastic/')
    try:
        # Run the command and wait for it to complete
        result = subprocess.run(command, check=True)

        # Check if the command was successful
        if result.returncode == 0:
            print("   ------>   Elastic build script executed successfully.")
        else:
            print(f"Elastic build script failed with return code: {result.returncode}")

    except subprocess.CalledProcessError as e:
        # Handle errors in execution
        print(f"Error running the elastic build script: {e}")

    os.chdir("../")
    
def run_elastic_clean(): 

    os.chdir('./elastic')
    command = ['./clean']  # Ensure there are no empty strings in the list

    try:
        # Run the command and wait for it to complete
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Print command output
        print(result.stdout)
        if result.returncode == 0:
            print("   ------>   Elastic clean script executed successfully.")
        else:
            print(f"Elastic clean script failed with return code: {result.returncode}")

    except FileNotFoundError:
        print("Error: 'clean' script not found in the 'elastic' directory.")
    except PermissionError:
        print("Error: 'clean' script is not executable. Try 'chmod +x ./clean'.")
    except subprocess.CalledProcessError as e:
        # Handle errors in execution
        print(f"Error running the elastic clean script: {e}")
        print(f"Standard Error:\n{e.stderr}")

    os.chdir("../")
    
     
def run_elastic_extract():
    """
    Function to run the command `python Elastic.git/Src/elastic.py extract`
    using subprocess.run.
    """
    command = ['python', 'Elastic.git/Src/elastic.py', 'extract']
    os.chdir('./elastic/')
    try:
        # Run the command and wait for it to complete
        result = subprocess.run(command, check=True)

        # Check if the command was successful
        if result.returncode == 0:
            print("   ------>   Elastic extract script executed successfully.")
        else:
            print(f"Elastic extract script failed with return code: {result.returncode}")

    except subprocess.CalledProcessError as e:
        # Handle errors in execution
        print(f"Error running the elastic extract script: {e}")

    os.chdir("../")
    
    
import os
import subprocess

# Function to detect all directories (bulk, babulk, c11, c44, etc.) and their subdirectories (10000, 10001, etc.)
def detect_directories_and_run_commands(tester_data, get_a0_data, extra_pot_files=None, elastic_dir="./elastic"):
    """
    Detect directories like bulk, babulk, c11, c44, and their subdirectories (10000, 10001, etc.),
    and execute the command `mpirun -np 4 EXE < cube.in > out.run` in each directory or subdirectory.

    - tester_data: Dictionary containing values from tester.input.
    - get_a0_data: Dictionary containing values from get_a0.input.
    - extra_pot_files: List of extra potential files to symlink.
    - elastic_dir: Path to the elastic directory.
    """
    # Get mpi_command and lammps_executable
    mpi_command = get_a0_data.get("mpi_command").split('#')[0].strip()  # Get mpi command from get_a0.input
    lammps_executable = tester_data.get("lammps_executable")  # Get lammps_executable from tester.input

    # Directories to search for
    main_dirs = ["bulk", "babulk", "c11", "c44"]

    for main_dir in main_dirs:
        main_path = os.path.join(elastic_dir, main_dir)

        # Check if the main directory exists
        if os.path.exists(main_path):
            # Check if this directory directly contains a cube.in file
            cube_in_file = os.path.join(main_path, "cube.in")
            if os.path.exists(cube_in_file):
                # Run the command in the main directory directly
                run_command_in_directory(main_path, mpi_command, lammps_executable, extra_pot_files)
                os.chdir("../..")
            else:
                # If no cube.in is found in the main directory, check for subdirectories
                subdirs = [d for d in os.listdir(main_path) if os.path.isdir(os.path.join(main_path, d))]

                for subdir in subdirs:
                    subdir_path = os.path.join(main_path, subdir)

                    # Check if cube.in exists in the subdirectory
                    cube_in_file = os.path.join(subdir_path, "cube.in")
                    if not os.path.exists(cube_in_file):
                        print(f"Warning: cube.in not found in {subdir_path}. Skipping.")
                        continue

                    # Run the command in each subdirectory
                    run_command_in_directory(subdir_path, mpi_command, lammps_executable, extra_pot_files)
                    os.chdir("../../..")


def run_command_in_directory(directory, mpi_command, lammps_executable, extra_pot_files=None):
    """
    Function to run the command in a specific directory.
    - directory: The directory where the command should be run.
    - mpi_command: The MPI command.
    - lammps_executable: The LAMMPS executable.
    - extra_pot_files: List of extra potential files to symlink.
    """
    command = mpi_command.split() + [lammps_executable]

    # Change to the specified directory and run the command
    os.chdir(directory)

    # Symlink extra potential files if provided
    if extra_pot_files:
        for f in extra_pot_files:
            if not os.path.exists(f):
                print(f"WARNING: extra potential file '{f}' not found, skipping")
                continue
            dest = os.path.basename(f)
            if not os.path.exists(dest):
                os.symlink(f, dest)

    try:
        #print(f"Executing in {directory}: {' '.join(command)}")

        out_run = "out.run"
        cube_in = "cube.in"
        #print(command, os.getcwd())
        with open(cube_in, 'r') as cube_in_py, open(out_run, 'w') as out_run_py:
            _ = subprocess.run(command,  stdin=cube_in_py, stdout=out_run_py, stderr=subprocess.PIPE)

        # print(f"   ------>   Command executed successfully in {directory}.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command in {directory}: {e}")

            
# Main function
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser('GetA0')
    parser.add_argument('-m','--mode',default="a0")
    args = parser.parse_args()
    mode = args.mode
    
    
    if mode == 'a0': 
      # Paths to input files
      tester_input_path = "../tester.input"  # Path to tester.input
      get_a0_input_path = "./get_a0.input"  # Path to get_a0.input
      # Output path for the LAMMPS input file
      output_lammps_input_path = "./elastic/LAMMPS/lammps_input.in"
      # Read input files
      tester_data = read_tester_input(tester_input_path)
      get_a0_data = read_get_a0_input(get_a0_input_path)
      # Validate inputs
      if not validate_inputs(tester_data, get_a0_data):
          sys.exit(1)
      
      # Extra potential files (optional)
      extra_pot_files = []
      for key in ['potential_file_ext1', 'potential_file_ext2']:
          val = tester_data.get(key, "").strip()
          if val:
              extra_pot_files.append(val)
      if extra_pot_files:
          print(f"   ------>   extra_pot_files   = {extra_pot_files}")
      
      
      # Create the LAMMPS input file
      create_lammps_input_file(tester_data, output_lammps_input_path)
      # Modify elastic/input_elastic_ini_org and save to elastic/input_elastic_ini
      input_elastic_org_path = "./elastic/input_elastic.ini_org"
      output_elastic_ini_path = "./elastic/input_elastic.ini"
      modify_elastic_input(get_a0_data, tester_data, input_elastic_org_path, output_elastic_ini_path)
      # Run the Elastic build script
      run_elastic_build()
      
      # Detect directories and run commands
      detect_directories_and_run_commands(tester_data, get_a0_data, extra_pot_files, elastic_dir="./elastic")
      
      # Run the Elastic extract script
      run_elastic_extract()
    elif mode == "clean": 
      run_elastic_clean()  
    else :
      raise NotImplementedError(f'Mode : {mode} is not implemented')
    

