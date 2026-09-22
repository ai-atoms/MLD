&input_ml
debug=.false.

ml_type=-1              !set -1 to compute descriptors only
write_desc=.true.       !set true to write the files to descDB/
desc_forces=.false.     !set true to compute the descriptors of forces


!Define your system
weighted=.false.           !set true for multicomponent systems
chemical_elements=" W "   !provide the system composition


!Descriptor settings
r_cut=4.7d0                !set the cutoff distace Rc
descriptor_type=4          !set 4 for bispectrum AFS
afs_type=1
n_rbf=11
n_cheb=11
&end
