&input_ml
debug=.false.

ml_type=-1              !set -1 to compute descriptors only
write_desc=.true.       !set true to write the files to descDB/
desc_forces=.false.     !set true to compute the descriptors of forces


!Define your system

weighted=.true. 
weighted_auto=.true.
weighted_3ch=.true.
fix_no_of_elements=4
chemical_elements=" Ta Ti V W  "
weight_per_element="0.1E+01  0.1E+01  0.1E+01  0.1E+01 "
weight_per_element_3ch="0.1E+01  0.1E+01  0.1E+01  0.1E+01 "

!Descriptor settings
r_cut=4.7d0                !set the cutoff distace Rc
descriptor_type=9          !set 4 for bispectrum AFS
j_max=4.0
&end
