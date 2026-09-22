&input_ml
debug=.false.

ml_type=-1              !set -1 to compute descriptors only
write_desc=.true.       !set true to write the files to descDB/
desc_forces=.false.     !set true to compute the descriptors of forces


!Define your system
weighted=.false.           !set true for multicomponent systems
chemical_elements=" Fe "   !provide the system composition


!Descriptor settings
r_cut=5.0d0                !set the cutoff distace Rc
descriptor_type=9          !set 9 for bispectrum SO4
j_max=4.0                  !angular moment for bispectrum SO4
Nfix=.true.
delta_fix_N_rcut=0.8
discrete_fix_N_rcut=100
&end
