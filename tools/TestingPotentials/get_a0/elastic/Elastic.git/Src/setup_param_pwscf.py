import globalv
import os

def init():
 global PWX, PHX, Q2RX, MATDYNX, ELASTIC, SAUVE2GIN, eampot
 global incar, kpoints
 global dinfile_cg, dinfile_cg2,dinfile_tr
 global base, inpdir
 global mdsub_o, inpmab, namein, ginin 
 global ph_in, phrecover_in, q2r_in, matdyn_in
 global flmodes, fldos,eigenvalues_phondy
 global prefix, relaxation
 global limit_omega
 global ntasks, nodes, npools, nimages
 nrepbox=4.0 
 base=os.getcwd()
 inpdir=base+'/INP_NDM/' 


 root_dir=globalv.root_dir
 dirPWSCF=root_dir+'/PWSCF'
 #potcar=dirPWSCF + '/' + 'POTCAR.UPF'
 incar=dirPWSCF + '/' + 'INCAR.pwscf'
 kpoints=dirPWSCF + '/' + 'KPOINTS.pwscf'
 #prefix='f2bcc'
 prefix='febcc'
 #prefix='fec15'

# curie zelda occigen
 computer=globalv.cluster
 if computer=='occigen':
     nodes=51

 #ntasks=12
 #npools=6
 #nimage=2


 #bcc
 #ntasks=640 
 #npools=16
 #nimage=10



 #c15
 ntasks=800 
 npools=20


 # occigen, zelda, my_computer
 #zelda 
 if computer=='zelda':
    PWX='/soft/QE/QE-5.4.0/espresso-5.4.0/bin/pw.x'
 #curie
 if computer=='curie':
    PWX='/ccc/products/espresso-5.1.1/default/bin/pw.x'
 if computer=='occigen':
    PWX='/opt/software/applications/espresso/5.1.1/bin/pw.x'



 for filetest in (PWX,incar,kpoints):
    if not os.path.exists(filetest):
       print("setup_param_pwscf: file %s doesn't exist. Put the correct path in setup_param_pwscf"%str(filetest))
       exit(0)

 script=globalv.script
 cal_elas=globalv.cal_elas
 if cal_elas=='qha_bulk':
    #limit for zero in frequencies (THz untis).  In order to detect  the three PBC modes
    limit_omega=0.03
    relaxation='no'
    # -ni for phonons calculation ; or for neb
    #nimage=2
    #bcc 
    nimages=10
    #c15  nimages=10
    ph_in       =dirPWSCF + '/' + 'PH.pwscf'
    phrecover_in=dirPWSCF + '/' + 'PHrecover.pwscf'
    q2r_in      =dirPWSCF + '/' + 'Q2R.pwscf'
    matdyn_in   =dirPWSCF + '/' + 'MATDYN.pwscf'
    #should be nodes multiply by 24 ntasks-per-node
    eigenvalues_phondy="eigenvalues.dat"
    if computer=='zelda':
      PHX='/soft/QE/QE-5.4.0/espresso-5.4.0/bin/ph.x'
      Q2RX='/soft/QE/QE-5.4.0/espresso-5.4.0/bin/q2r.x'
      MATDYNX='/soft/QE/QE-5.4.0/espresso-5.4.0/bin/matdyn.x'
    if computer=='curie':
      PHX='/ccc/products/espresso-5.1.1/default/bin/ph.x'
      Q2RX='/ccc/products/espresso-5.1.1/default/bin/q2r.x'
      MATDYNX='/ccc/work/cont002/den/marinica/QE-5.4.0/espresso-5.4.0/bin_curie/matdyn.x'
    if computer=='occigen':
      PHX='/opt/software/applications/espresso/5.1.1/bin/ph.x'
      Q2RX='/opt/software/applications/espresso/5.1.1/bin/q2r.x'
      MATDYNX='/store/mcm/QE-5.4.0/espresso-5.4.0/bin/matdyn.x'

    for testf in (PHX, Q2RX, MATDYNX, ph_in, phrecover_in, q2r_in, matdyn_in):
      if not os.path.exists(testf):
        print("ndm setup_ndm_pwscf: The file %s doesn't exist. Put the correct path. "%(str(testf)))
        exit(0)




 return ;
