
############################################# Required Modules and Libraries##########################

import numpy as np
import pandas as pd
import glob
from Compute_delta_properties import Compute_change_in_time_Grain, Compute_change_in_time_GrainBoundary

from MicroTANA_class import Microstructure


###################### Input Files  ####################################

dir_path="E:/Anneal500+Compression50Al161/GradeAcalcA500C50Al161"#"/home/uda69/lcz56/SimulationDataFromOccigen/Al/Al_161x161x161/LAMMPS/annealing500Al161x161x161/Quenched_After_AnnealAl161_500P500/GradeA_Al161_500P500"#"/home/uda69/lcz56/GradeA_Al250_500"
core_name ="dump.anneal500Al161_P50"# "dump.annealAl161P500"
nf=97
n_analyze=2
filename=[]
for i in range (1,nf):
    
    a=dir_path+"/"+core_name+"__AtomData_" + str(i)+ "_*.cfg"
    print(a)
    filename.append(glob.glob(a)[0])

file_path=dir_path+"/"+core_name+"_grain_id_mapping.npy"

grain_mapping_data =np.load(file_path)

# The following are file paths for  specifying the orientation of each grains at different time-steps 
# The files are outputes of Grade-A

file_path_q0=dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_q0.csv"
file_path_q1=dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_q1.csv"
file_path_q2=dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_q2.csv"
file_path_q3=dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_q3.csv"
file_path_cubMis=dir_path+"/TimeEvo/"+core_name+"__GrainData__TimeEvo_cubMisOri.csv"


q0=pd.read_csv(file_path_q0,sep=';')
q0=(pd.DataFrame(q0).to_numpy())[:,1:]
q1=pd.read_csv(file_path_q1,sep=';')
q1=(pd.DataFrame(q1).to_numpy())[:,1:]

q2=pd.read_csv(file_path_q2,sep=';')
q2=(pd.DataFrame(q2).to_numpy())[:,1:]
q3=pd.read_csv(file_path_q3,sep=';')
q3=(pd.DataFrame(q3).to_numpy())[:,1:]
cubMis=pd.read_csv(file_path_cubMis,sep=';')
cubMis=(pd.DataFrame(cubMis).to_numpy())[:,1:]


################### Helper Functions #################################
def get_Micro_ordered_topology(M_info):
    
    filename=M_info['filename']
    timestep=M_info['timestep']
    ThresholdDistance=M_info['ThresholdDistance']
    
    isPeriodic=M_info['isPeriodic']
    atomic_radius=M_info['atomic_radius']
    
    M=Microstructure(filename,timestep,ThresholdDistance,isPeriodic,atomic_radius)
    M.findAtomsGrain_c()
    M.PeriodicNeighbourIndex_c()
    M.UpdateGrainIdForUnallocatedOrphanAtoms_c()
    M.ModifiedGrainIdComputation_c()
    M.Update_x_c()
    M.Ordered_Topologies_c()
    
    return M
    
def Run_Micro(M,G_lst,grain_mapping_data,cubMis,method):
    M.construct_triplelines(triplelineID_list=[])

    #G_lst=[(0,),(1,),(14,)]
    M.construct_grains(method=method,grainID_list=G_lst)



    GB_lst=[]

    for G_i in G_lst:
        
        GB_lst=GB_lst+[GB_id for GB_id in M.grainboundary_id_list if set(G_i).issubset(GB_id)]

    GB_lst_temp=[GB_i for GB_i in GB_lst if  False==((GB_i[0],) in G_lst and (GB_i[1],) in G_lst)]

    G_lst2=[tuple([i for i in GB_i if ((i,) not in G_lst) ]) for GB_i in GB_lst_temp]

    M.construct_grains(method=method,grainID_list=G_lst2)

    M.construct_grainboundaries(method=method,tolerance=[],grainboundaryID_list=GB_lst)


    for GB_id,GB in M.grainboundaries.items():
        GB.compute_misorientation((q0,q1,q2,q3),grain_mapping_data)
        GB.compute_imc_surfacearea_volume_GB()
        GB.compute_topological_properties_GB(M=M)
     
    for G_id,G in M.grains.items():
        G.compute_imc_surfacearea_volume_G()
        G.compute_topological_properties_G(M=M)
        G.compute_relative_rotation(M,grain_mapping_data, cubMis)
    
    G_lst_combined=G_lst+G_lst2
    return (M,G_lst_combined,GB_lst)

    

#############  Analysis  Microstructure #######################################

t=[0,1]
M_info0={'filename':filename[0],'timestep':t[0],'ThresholdDistance':3.25,'isPeriodic':True,
         'atomic_radius':1.43}
M0=get_Micro_ordered_topology(M_info0)

M_info1={'filename':filename[1],'timestep':t[1],'ThresholdDistance':3.25,'isPeriodic':True,
         'atomic_radius':1.43}
M1=get_Micro_ordered_topology(M_info0)

G_lst=[(0,),(1,),(14,)]
M0,G_lst_combined,GB_lst=Run_Micro(M=M0,G_lst=G_lst,grain_mapping_data=grain_mapping_data,cubMis=cubMis,method=1)

Output_G=[]
for G_id0 in G_lst:

    Out_G=Compute_change_in_time_Grain(M0,M1,G_id0,grain_mapping_data,cubMis,Method='ray')
    Output_G.append(Out_G)

Output_GB=[]    
for GB_id0 in GB_lst:
    Out_GB=Compute_change_in_time_GrainBoundary(M0,M1,GB_id0,grain_mapping_data,Use_grain=True)
    Output_GB.append(Out_GB)
    
  
    



    
