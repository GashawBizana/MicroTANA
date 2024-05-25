# MicroTANA
MicroTANA is a python program that can be used to extract information that are relevant to grain growth from Molecular dynamics trajectory data.
The program relies on the following libraries and programs:Numpy, Scipy, Pandas, Ovito, trimesh, pyvisat,pymeshfix,  embree, and Grade-A (https://github.com/paulhof/GraDe-A)
 The program is run using the Microstructure Class  in MicroTANA_class.py. A script showing usage of the program given in the script MT_cb.py . 
 Note that the program requires the outputs of the CFG files and Time-Evolution files that are obatined from Grade-A and grain mapping file that is created using the outputs of Grade-A.
 THE grain_mapping_data can be created from the outputs of Grade-A using the GrainId_mapping.py script or the matlab code GrainIdCreat.m.
