import os
import re
import sys
import numpy as np
import subprocess 

def run_ndm(mld_exe, nprocs):

  #os.system("mpirun -np %s %s -> ndm.out" % (nprocs, mld_exe))

  command = ["mpirun", "-np", "%s"%nprocs, "%s"%mld_exe]
  result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

  # Check if the execution was successful
  if result.returncode == 0:
    #print("Execution successful.")
    #print("Output:", result.stdout.decode())
    with open('ndm.out', 'w') as file:
        file.write(result.stdout.decode())

  else:
    print("Execution failed.")
    print("Error:", result.stderr.decode())

  return



def get_natom_from_ndm(outNDM):
 f = open(outNDM, 'r')
 word = []
 for line in f:
   if re.match(' nombre d atomes ', line):
     word = line.split()
     #print word[4]
 f.close()

 if len(word) == 0:
   print('get_natom_from_ndm: There is no nat output in the directory %s' % os.getcwd())
   exit(0)

 return float(word[4])


def get_energy_from_ndm(outNDM):
 f = open(outNDM, 'r')
 word = []
 for line in f:
   if re.match('(.*) energie (.*)', line):
     word = line.split()
     #print word
 f.close()

 if len(word) == 0:
   print('get_energy_from_ndm: There is no energie output in the directory %s' % os.getcwd())
   exit(0)
 return float(word[-1])


def compare_energies(ref_file, cur_file):
   ref_energy = float(get_energy_from_ndm(ref_file))
   cur_energy = float(get_energy_from_ndm(cur_file))
   diff = abs(ref_energy - cur_energy)
   ltest = True
   if (diff > 1.e-10):
       ltest = False
   return ltest


def read_x(filename, i):
    x = np.loadtxt(filename, usecols=(i))

    return x



def compare_fit(ref_DIR, cur_DIR):
   ref_ene = read_x(ref_DIR + '/' + 'train_energy.out')
   fit_ene = read_x(cur_DIR + '/' + 'train_energy.out')
   diff = abs(ref_energy - cur_energy)
   ltest = True
   if (diff > 1.e-10):
       ltest = False
   return ltest

def check_path(path):
    if os.path.isdir(path) and not os.path.islink(path):
        return 1
    elif os.path.islink(path):
        return 2
    else:
        return 0

def process_files(directory_path, desired_columns):
    # 1. List all files in the directory
    statuss=check_path(directory_path)
    if (statuss == 0):
      print(directory_path + ' does not exist')
    all_files = os.listdir(directory_path)
    
    eml_files = [file for file in all_files if file.endswith('.eml')]
    
    # Store content based on number of lines
    content_dict = {}

    for file in eml_files:
        file_path = os.path.join(directory_path, file)   
        # 2. Open each file and check columns
        with open(file_path, 'r') as f:
            content = [line.split() for line in f]
            
            if len(content[0]) == desired_columns:
                num_lines = len(content)
            
                
                if num_lines not in content_dict:
                    content_dict[num_lines] = []
                content_dict[num_lines].append(content)
            else: 
                print(f'File {file} has {len(content[0])} columns instead of {desired_columns}')
                continue    
    # 3 & 4. Processing each group of files with same number of lines

    results = {}
    for num_lines, matrices in content_dict.items():
        concatenated = np.vstack(matrices)
        
        norms = np.linalg.norm(concatenated[:, 1:], axis=1)
        average_norm = np.mean(norms)
        
        results[num_lines] = average_norm
        
    # 5. Sort and convert to 2D numpy array
    sorted_results = sorted(results.items(), key=lambda x: x[0])
    results_array = np.array(sorted_results)
    
    return results_array 

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

#   energy=property(fget=_get_energy, fset=_set_energy, fdel=None, doc=None)

def compare_fit(fit_limit) :
  
  cur_fit = Fit()
  ref_fit = Fit()
  cur_fit.energy_train = 'train_energy.out'
  cur_fit.energy_dft = 'train_energy.out'
  ref_fit.energy_train = 'REF/train_energy.out'
  ref_fit.energy_dft = 'REF/train_energy.out'
  corr = np.corrcoef(cur_fit.energy_train, ref_fit.energy_train)[1, 0]
  print('                 .... correlation {:1.5f}'.format(corr))
  if (abs(corr - 1.0) < fit_limit):
     print('                 .... test energy OK')
  else:
     print('                 .... test energy FAILED')  
  cur_fit.force_train = 'train_force.out'
  cur_fit.force_dft = 'train_force.out'
  ref_fit.force_train = 'REF/train_force.out'
  ref_fit.force_dft = 'REF/train_force.out'
  corr = np.corrcoef(cur_fit.force_train, ref_fit.force_train)[1, 0]
  print('                 .... correlation {:1.5f}'.format(corr))
  if (abs(corr - 1.0) < fit_limit):
     print('                 .... test force OK')
  else:
     print('                 .... test force FAILED')  
  cur_fit.stress_train = 'train_stress.out'
  cur_fit.stress_dft = 'train_stress.out'
  ref_fit.stress_train = 'REF/train_stress.out'
  ref_fit.stress_dft = 'REF/train_stress.out'
  corr = np.corrcoef(cur_fit.stress_train, ref_fit.stress_train)[1, 0]
  print('                 .... correlation {:1.5f}'.format(corr))
  if (abs(corr - 1.0) < fit_limit):
     print('                 .... test stress OK')
  else:
     print('                 .... test stress FAILED')  

  return ;

def compare_write_desc(vout, vref):

  for i in range(len(vout)):
    ff=1e-10
    vdiff=abs(vout[i][1] - vref[i][1])
    if (abs(vdiff) < ff):
      print(f"               .... diff {int(vout[i][0]):>5}   OK")
    else:
      print(f"               .... diff {int(vout[i][0]):>5}   FAILED  {vdiff:.18f}")

  return ; 
#HERE THE SCRIPT ...................
usage = """

           Usage: python small_test.py mode

           mode    -  'run'  do the  tests
                      'clean'   clean the directories  ;


"""


if (len(sys.argv) != 2):
   print(usage)
   exit(0)

lrun = False
lclean = False
if sys.argv[1] == 'run':
    lrun = True

if sys.argv[1] == 'clean':
    lclean = True

if (not lclean) and (not lrun):
   print("The options are run or clean")
   print(usage)
   exit(0)

#mld_exe = os.environ['HOME'] + '/NDM/p2.8/rundm90_ml_para'
mld_exe='/volatile/home/az268520/GitHub/MLD/mld_build_gnu/bin/milady_main.exe'
#mld_exe = '/home/lapointe/Git/mld_dev_build_intel/bin/milady_main.exe'
nprocs = 2
BaseDir = os.getcwd()
fit_limit = 1.e-3
##Here are MD tests ....
"""
list_of_tests = ['NDM_2020/md_cg_ml',
                 'NDM_2020/md_minimize_ml']
for i in range(len(list_of_tests)):
   print(('Example MD %s     .... %s' % (i, list_of_tests[i])))
   os.chdir(BaseDir + '/' + list_of_tests[i])
   if lrun:
      run_ndm(mld_exe, nprocs)
      ltest = compare_energies('ref_out', 'ndm.out')
      if ltest:
         print('                 .... test OK')
      else:
         print('                 .... test FAILED')
   if lclean:
      os.system('./clean.sh')
"""

"""
list_of_tests = ['train_tests/fit19',
                 'train_tests/fit18',
                 'train_tests/fit17'
                 ]


"""


list_of_tests = ['train_tests/fit01',
                 'train_tests/fit02',
                 'train_tests/fit03',
                 'train_tests/fit04',
                 'train_tests/fit05',
                 'train_tests/fit06',
                 'train_tests/fit07',
                 'train_tests/fit08',
                 'train_tests/fit09',
                 'train_tests/fit10',
                 'train_tests/fit11',
                 'train_tests/fit13',
                 'train_tests/fit14',
                 'train_tests/fit15',
                 'train_tests/fit16',
                 'train_tests/fit17',
                 'train_tests/fit18',
                 'train_tests/fit19',
                 'train_tests/fit20',
                 'train_tests/fit21',
                 'train_tests/fit22',
                 'train_tests/fit23',
                 'train_tests/fit24',
                 'train_tests/fit25']



"""
list_of_tests = ['train_tests/fit16',
                 'train_tests/fit17',
                 'train_tests/fit18',
                 'train_tests/fit19',
                 'train_tests/fit20']
"""

if lrun:
  os.chdir(BaseDir + '/' +'train_tests/data_for_tests')
  os.system('tar -zxvf DB_HEA.tgz > /dev/null')
  os.system('tar -zxvf DB_one.tgz > /dev/null')
  os.system('tar -zxvf DB_two.tgz > /dev/null')
if lclean: 
  os.chdir(BaseDir + '/' +'train_tests/data_for_tests')
  os.system('rm -rf DB_HEA')  
  os.system('rm -rf DB_one')  
  os.system('rm -rf DB_two')  
  

for i in range(len(list_of_tests)):
   print(('Example TRAIN  %s .... %s' % (i, list_of_tests[i])))
   os.chdir(BaseDir + '/' + list_of_tests[i])
   if lrun:
      f = open('README', 'r')
      for line in f:
         print(line)
      f.close()
     
      if (list_of_tests[i] == 'train_tests/fit12'):
        os.system("cp vacancy.ml_write_dump  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_read_dump  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        #os.system("cp -f *.out REF/")
        #os.system("cp -f *.milady REF/")
      elif (list_of_tests[i] == 'train_tests/fit13'):
        os.system("cp vacancy.ml_01  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_02  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_03  vacancy.ml")
        run_ndm(mld_exe, nprocs)  
        #os.system("cp -f *.out REF/")
        #os.system("cp -f *.milady REF/")
      elif (list_of_tests[i] == 'train_tests/fit14'):
        os.system("cp vacancy.ml_01  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_02  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_03  vacancy.ml")
        run_ndm(mld_exe, nprocs) 
        #os.system("cp -f *.out REF/")
        #os.system("cp -f *.milady REF/")
     
      elif (list_of_tests[i] == 'train_tests/fit15'):
        os.system("cp vacancy.ml_01  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_02  vacancy.ml")
        run_ndm(mld_exe, nprocs)
        os.system("cp vacancy.ml_03  vacancy.ml")
        run_ndm(mld_exe, nprocs)       
        #os.system("cp -f *.out REF/")
        #os.system("cp -f *.milady REF/")
      
                 
      else:  
        run_ndm(mld_exe, nprocs)

        #os.system("cp -f *.out REF/")
        #os.system("cp -f *.milady REF/")
      
      if (list_of_tests[i] == 'train_tests/fit17'):
        vout = process_files('descDB/', 56)
        vref = process_files('REF/descDB/', 56)
        #os.system("rm  -f REF/descDB/*.eml")
        #os.system("cp descDB/*.eml REF/descDB/")
        compare_write_desc(vout, vref)
      elif (list_of_tests[i] == 'train_tests/fit18'):
        vout = process_files('descDB/', 133)
        vref = process_files('REF/descDB/', 133)
        #os.system("rm  -f REF/descDB/*.eml")   
        #os.system("cp descDB/*.eml REF/descDB/")
        compare_write_desc(vout, vref)
      elif (list_of_tests[i] == 'train_tests/fit19'):
        vout = process_files('descDB/', 166)
        vref = process_files('REF/descDB/', 166)
        #os.system("rm  -f REF/descDB/*.eml")
        #os.system("cp descDB/*.eml REF/descDB/")
        compare_write_desc(vout, vref)

      else:
         compare_fit(fit_limit)

   if lclean:
      os.system('./clean.sh')

   #os.system('pwd')
