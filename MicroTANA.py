###### External Modules ######
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import sys
import random
import time
import concurrent.futures
import importlib.util
import glob
from multiprocessing import Pool
import multiprocessing
from contextlib import closing
import scipy.io
import gc
########## Import Local Modules ########
# import Data_IO 
# import UnallocatedOrphanAtomGrainIdUpdater 
# import PeriodicNearestNeighbourFinder 
# import ModifiedGrainIdIdentifier 
# import GrainDataConstructor 

def PeriodicDistance(X1,X21,EdgeLength,isPeriodic):
    d=[]
    for i in range(0,len(X21)):
        X2=X21[i]
        
        dx=X1[0]-X2[0]
        dy=X1[1]-X2[1]
        dz=X1[2]-X2[2]
        
        if (isPeriodic==1):
            if (dx>EdgeLength[0]*0.5):
                dx=dx-EdgeLength[0]
            elif(dx<=-EdgeLength[0]*0.5):
                dx=dx+EdgeLength[0]
        
        if (isPeriodic==1):
            if (dy>EdgeLength[1]*0.5):
                dy=dy-EdgeLength[1]
            elif(dy<=-EdgeLength[1]*0.5):
                dy=dy+EdgeLength[1]
        
        if (isPeriodic==1):
            if (dz>EdgeLength[2]*0.5):
                dz=dz-EdgeLength[2]
            elif(dz<=-EdgeLength[2]*0.5):
                dz=dz+EdgeLength[2]
        
        d.append(np.sqrt(np.square(dx)+np.square(dy)+np.square(dz)))
    return(d)

def PeriodicCenterOfMass(xn,EdgeLength):
    
    thetai=2*np.pi*xn/EdgeLength
    etai=np.cos(thetai)
    gammai=np.sin(thetai)
    etai_mean=np.mean(etai,axis=0)
    gammai_mean=np.mean(gammai,axis=0)
    theta_mean=np.arctan2(-gammai_mean,-etai_mean)+np.pi
    xcm=EdgeLength*(theta_mean/(2*np.pi))
    return(xcm.tolist())


##### Internal Modules #############
from Data_IO import GrainEdgeLength,ReadGrainData,OutPutCFGWritter, findAtomsGrain
from UnallocatedOrphanAtomGrainIdUpdater import UpdateGrainIdForUnallocatedOrphanAtoms
from PeriodicNearestNeighbourFinder import PeriodicNeighbourIndex
from ModifiedGrainIdIdentifier import ModifiedGrainIdComputation
from GrainDataConstructor import OrderedTopologies
from NormalEstimation import NormalvectorEstimationWhole
from Curvature_analysis import GrainBoundaryProperties,GrainBoundaryPropertyExtended
from Curvature_analysis_ovito import GrainBoundaryProperties_ovito,GrainBoundaryPropertyExtended_ovito
from TripleLineLength_computation import TripleLineLength_For_all
from Microstructure_Topology_info import Microstructure_topology
# from microstructure_mesh_and_property import GrainProperty,GrainBoundaryProperty
# from Velocity_computation import Compute_distance_and_average_curavture

####### Computation one files  ###########
def MicroTANA_one_file(file_and_timestep,ThresholdDistance=3.25,isPeriodic=True):
    ###### Input Data ##################
    file=file_and_timestep[0]
    timestep=file_and_timestep[1]
    X,GrainId,Image,H0,n_p,EdgeLength=findAtomsGrain(file)
    ImageShift=Image*EdgeLength
    
    print (file + " MicroTANA Starting ...")
    t = time.time()
    
    #####  Periodic Neighbour List Tree #########
    t1=time.time()
    print (file + " Importing Data took:", (t1-t))
    
    print (file + " Creating Neighbour List ...")
    PeriodicNeighbourList_Index=PeriodicNeighbourIndex(EdgeLength,X,ThresholdDistance)
    t2=time.time()
    print (file + " Creating Neighbour List took:", (t2-t1))
    ########## Update Grade-A loose atom  Grain Id #######
    
    
    print (file + " Adopting Loose atoms ...")
    
    ### Exaustively use the grain adoption of GradeA
    UpdateGrainIdForUnallocatedOrphanAtoms(PeriodicNeighbourList_Index,GrainId,X,EdgeLength,isPeriodic)
    t3=time.time()
    print (file + " Adopting Loose atoms took:" , (t3-t2))
    ###### Modified GrainId List #############
    print (file + " Assigning Modified Grain Id ...")
    
    ModifiedGrainIdListSorted=ModifiedGrainIdComputation(PeriodicNeighbourList_Index,GrainId,X,EdgeLength,isPeriodic)
    t4=time.time()
    print (file + " Assigning Modified Grain Id  took:", (t4-t3))
    #ModifiedGrainIdListSorted=sorted(ModifiedGrainIdList, key=lambda x: x[0]) ## use this if mpi parallizing
    
    Xnew=[]
    for i,x in enumerate(X):
        gid=GrainId[i]
        mgid=ModifiedGrainIdListSorted[i][1]
        tl=ModifiedGrainIdListSorted[i][2]
        PNL=PeriodicNeighbourList_Index[i]
        if mgid >=0:
            Xnew.append(x.tolist())
        elif mgid ==-1:
            gid2=[j for j in tl if j!=gid]
            xn1=[X[j] for j in PNL if ((ModifiedGrainIdListSorted[j][1]==-1 and GrainId[j] in gid2) or (X[j]==x).all())]
            dist=PeriodicDistance(x,xn1,EdgeLength,isPeriodic)
            idxi=sorted(range(len(dist)), key=lambda k: dist[k])[0:2]
            Xnew.append(PeriodicCenterOfMass(np.array(xn1)[idxi],EdgeLength))
            
        elif mgid ==-2:
            gid2=[j for j in tl if j!=gid]
            xn1=[X[j] for j in PNL if ((ModifiedGrainIdListSorted[j][1]==-2 and GrainId[j] in gid2) or (X[j]==x).all())]
            dist=PeriodicDistance(x,xn1,EdgeLength,isPeriodic)
            idxi=sorted(range(len(dist)), key=lambda k: dist[k])[0:3]
            Xnew.append(PeriodicCenterOfMass(np.array(xn1)[idxi],EdgeLength))
            
        elif mgid ==-3:
             gid2=[j for j in tl if j!=gid]
             xn1=[X[j] for j in PNL if ((ModifiedGrainIdListSorted[j][1]==-3 and GrainId[j] in gid2) or (X[j]==x).all())]
             dist=PeriodicDistance(x,xn1,EdgeLength,isPeriodic)
             idxi=sorted(range(len(dist)), key=lambda k: dist[k])[0:4]
             Xnew.append(PeriodicCenterOfMass(np.array(xn1)[idxi],EdgeLength))
        else:
            Xnew.append(x.tolist())
    Xnew=np.array(Xnew)

    ###### Structured Data for Microstructure Newwork entities #####
    
    print (file + " Creating Microstructure Network Entities ...")
    InnerAtom,GrainBoundary,TripleLine,QuadraplePoint,Curvature=OrderedTopologies(ModifiedGrainIdListSorted,Xnew)
    t5=time.time()
    print (file + " Creating Microstructure Network Entities took: ",(t5-t4))
    ####### Write Modified Grain Id Output into CFG file ########
    print (file + " Writting Output files with Modified Grain Id ...")
    
    # OutPutCFGWritter(ModifiedGrainIdListSorted,GrainData,GrainId,file,fname)
    t6=time.time()
    print (file + " Writting Output files with Modified Grain Id took: ", (t6-t5))
    ###################### Normal Estimation ################
    
    # NormalVector=NormalvectorEstimationWhole(PeriodicNeighbourList_Index,X,EdgeLength)
    print (file + " Curvature analysis  ...")
   
    grain_boundary_property_extended_current,grain_property_current=GrainBoundaryProperties_ovito(GrainBoundary,InnerAtom,EdgeLength,Xnew,GrainId,ModifiedGrainIdListSorted,TripleLine,QuadraplePoint,atomicradii=1.43,method=1)
    grain_boundary_property_current=GrainBoundaryPropertyExtended_ovito(GrainBoundary,grain_boundary_property_extended_current)
    triple_line_length_current=TripleLineLength_For_all(TripleLine,EdgeLength)
    Grain_number_of_grain_boundaries_current, Grain_number_of_triple_lines_current,Grain_average_misorientation_current,Grain_length_of_triple_lines_current,grain_total_curvature_surface_area_volume_current, GrainBoundary_number_of_triple_line_current,GrainBoundary_misorientation_current,GrainBoundary_length_of_triple_line_current,GrainBoundary_total_curvature_mean_curvature_surface_area_current,grain_total_curvature_surface_area_volume_from_face_current=Microstructure_topology(InnerAtom,GrainBoundary,TripleLine,grain_property_current,grain_boundary_property_current, grain_boundary_property_extended_current,triple_line_length_current,timestep,q0,q1,q2,q3,grain_mapping_data,Xnew,GrainId,EdgeLength)
    
    t7=time.time()
    print (file + " Curvature analysis for one time step took took: ", (t7-t6))
    
    print (file + " MicroTANA took: ", (time.time()-t))
    
    return([grain_boundary_property_current, grain_property_current,triple_line_length_current,Grain_number_of_grain_boundaries_current, Grain_number_of_triple_lines_current,Grain_average_misorientation_current,Grain_length_of_triple_lines_current,grain_total_curvature_surface_area_volume_current, GrainBoundary_number_of_triple_line_current,GrainBoundary_misorientation_current,GrainBoundary_length_of_triple_line_current,GrainBoundary_total_curvature_mean_curvature_surface_area_current,grain_total_curvature_surface_area_volume_from_face_current,grain_boundary_property_extended_current, EdgeLength,InnerAtom,GrainBoundary,TripleLine,QuadraplePoint])

''' out put og microTANA is
0- grain_boundary_property, 1- grain_property, 2- triple_line_length, 3- Grain_number_of_grain_boundaries, 4-Grain_number_of_triple_lines, 5- Grain_average_misorientation

6-Grain_length_of_triple_lines, 7- grain_total_curvature_surface_area_volume, 8 - GrainBoundary_number_of_triple, 9- GrainBoundary_misorientation, 10 -GrainBoundary_length_of_triple_line,

11- GrainBoundary_total_curvature_surface_area_volume_current, 11- EdgeLength of simulation box

'''

####### Computation for multiple files  ###########

####### file names   ###########

filename=[]
for i in range (1,65):
   
    a="F:/GradeA_Al161_500/dump.annealAl161__AtomData_" + str(i)+ "_*.cfg"
    
    filename.append(glob.glob(a)[0])
file_path="F:/GradeA_Al161_500/A500Al161_grain_id_mapping.mat"
mat = scipy.io.loadmat(file_path)
grain_mapping_data = mat['TrackedId']

file_path_q0='F:/GradeA_Al161_500/TimeEvo/dump.annealAl161__GrainData__TimeEvo_q0.csv'
file_path_q1='F:/GradeA_Al161_500/TimeEvo/dump.annealAl161__GrainData__TimeEvo_q1.csv'
file_path_q2='F:/GradeA_Al161_500/TimeEvo/dump.annealAl161__GrainData__TimeEvo_q2.csv'
file_path_q3='F:/GradeA_Al161_500/TimeEvo/dump.annealAl161__GrainData__TimeEvo_q3.csv'

q0=pd.read_csv(file_path_q0,sep=';')
q0=(pd.DataFrame(q0).to_numpy())[:,1:]
q1=pd.read_csv(file_path_q1,sep=';')
q1=(pd.DataFrame(q1).to_numpy())[:,1:]

q2=pd.read_csv(file_path_q2,sep=';')
q2=(pd.DataFrame(q2).to_numpy())[:,1:]
q3=pd.read_csv(file_path_q3,sep=';')
q3=(pd.DataFrame(q3).to_numpy())[:,1:]

######## Multi processing #########

if __name__ == '__main__':
    ii=[ii for ii in range(0,3)]
    a=filename[0:3]
    b= [i for i in range(0,3)]
    aa=[]
    for i in range(0,len(a)):
        aa.append([a[i],ii[i]])
    
    npr=int(0.5*multiprocessing.cpu_count())
    cs=len(aa)/npr
    
    grain_boundary_property_for_all_timesteps, grain_property_for_all_timesteps,triple_line_length_for_all_timesteps,Grain_number_of_grain_boundaries_for_all_timesteps, Grain_number_of_triple_lines_for_all_timesteps,Grain_average_misorientation_for_all_timesteps,Grain_length_of_triple_lines_for_all_timesteps,grain_total_curvature_surface_area_volume_for_all_timesteps, GrainBoundary_number_of_triple_line_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps,GrainBoundary_length_of_triple_line_for_all_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps,grain_total_curvature_surface_area_volume_from_face_for_all_timesteps,grain_boundary_property_extended_for_all_timesteps, EdgeLength_for_all_timesteps,InnerAtom_for_all_timesteps,GrainBoundary_for_all_timesteps,TripleLine_for_all_timesteps,QuadraplePoint_for_all_timesteps=(),(),(),(),(),(),(),(),(),(),(),(),(),(),(),(),(),(),()
    
    for idx in range(0,math.ceil(cs)):
        if idx==0:
            inputfiles=aa[idx*npr:(idx+1)*npr]
            # print(inputfiles) 
            
        else:
            inputfiles=aa[idx*npr:(idx+1)*npr]
            # print(inputfiles)
        
        if not inputfiles:
            break
        
        a_pool =Pool(npr)
    

        ''' out put og microTANA is
        0- grain_boundary_property, 1- grain_property, 2- triple_line_length, 3- Grain_number_of_grain_boundaries, 4-Grain_number_of_triple_lines, 5- Grain_average_misorientation
    
        6-Grain_length_of_triple_lines, 7- grain_total_curvature_surface_area_volume, 8 - GrainBoundary_number_of_triple, 9- GrainBoundary_misorientation, 10 -GrainBoundary_length_of_triple_line,
    
        11- EdgeLength
    
        '''
    
        grain_boundary_property_for_chunk_timesteps, grain_property_for_chunk_timesteps,triple_line_length_for_chunk_timesteps,Grain_number_of_grain_boundaries_for_chunk_timesteps, Grain_number_of_triple_lines_for_chunk_timesteps,Grain_average_misorientation_for_chunk_timesteps,Grain_length_of_triple_lines_for_chunk_timesteps,grain_total_curvature_surface_area_volume_for_chunk_timesteps, GrainBoundary_number_of_triple_line_for_chunk_timesteps,GrainBoundary_misorientation_for_chunk_timesteps,GrainBoundary_length_of_triple_line_for_chunk_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_chunk_timesteps,grain_total_curvature_surface_area_volume_from_face_for_chunk_timesteps,grain_boundary_property_extended_for_chunk_timesteps, EdgeLength_for_chunk_timesteps,InnerAtom_for_chunk_timesteps,GrainBoundary_for_chunk_timesteps,TripleLine_for_chunk_timesteps,QuadraplePoint_for_chunk_timesteps = zip(*a_pool.imap(MicroTANA_one_file,inputfiles))
        ''' Need TO FIX MISORIENTATION FOR GB'''
        a_pool.close()
        a_pool.join()
        
        grain_boundary_property_for_all_timesteps, grain_property_for_all_timesteps,triple_line_length_for_all_timesteps,Grain_number_of_grain_boundaries_for_all_timesteps, Grain_number_of_triple_lines_for_all_timesteps,Grain_average_misorientation_for_all_timesteps,Grain_length_of_triple_lines_for_all_timesteps,grain_total_curvature_surface_area_volume_for_all_timesteps, GrainBoundary_number_of_triple_line_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps,GrainBoundary_length_of_triple_line_for_all_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps,grain_total_curvature_surface_area_volume_from_face_for_all_timesteps,grain_boundary_property_extended_for_all_timesteps, EdgeLength_for_all_timesteps,InnerAtom_for_all_timesteps,GrainBoundary_for_all_timesteps,TripleLine_for_all_timesteps,QuadraplePoint_for_all_timesteps=grain_boundary_property_for_all_timesteps + grain_boundary_property_for_chunk_timesteps,grain_property_for_all_timesteps + grain_property_for_chunk_timesteps,triple_line_length_for_all_timesteps + triple_line_length_for_chunk_timesteps,Grain_number_of_grain_boundaries_for_all_timesteps + Grain_number_of_grain_boundaries_for_chunk_timesteps,Grain_number_of_triple_lines_for_all_timesteps + Grain_number_of_triple_lines_for_chunk_timesteps,Grain_average_misorientation_for_all_timesteps + Grain_average_misorientation_for_chunk_timesteps,Grain_length_of_triple_lines_for_all_timesteps + Grain_length_of_triple_lines_for_chunk_timesteps,grain_total_curvature_surface_area_volume_for_all_timesteps + grain_total_curvature_surface_area_volume_for_chunk_timesteps,GrainBoundary_number_of_triple_line_for_all_timesteps + GrainBoundary_number_of_triple_line_for_chunk_timesteps,GrainBoundary_misorientation_for_all_timesteps + GrainBoundary_misorientation_for_chunk_timesteps,GrainBoundary_length_of_triple_line_for_all_timesteps + GrainBoundary_length_of_triple_line_for_chunk_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps + GrainBoundary_total_curvature_mean_curvature_surface_area_for_chunk_timesteps,grain_total_curvature_surface_area_volume_from_face_for_all_timesteps + grain_total_curvature_surface_area_volume_from_face_for_chunk_timesteps,grain_boundary_property_extended_for_all_timesteps + grain_boundary_property_extended_for_chunk_timesteps,EdgeLength_for_all_timesteps + EdgeLength_for_chunk_timesteps,InnerAtom_for_all_timesteps + InnerAtom_for_chunk_timesteps,GrainBoundary_for_all_timesteps + GrainBoundary_for_chunk_timesteps,TripleLine_for_all_timesteps + TripleLine_for_chunk_timesteps,QuadraplePoint_for_all_timesteps +QuadraplePoint_for_chunk_timesteps
        
        # print(EdgeLength_for_chunk_timesteps)

        del grain_boundary_property_for_chunk_timesteps, grain_property_for_chunk_timesteps,triple_line_length_for_chunk_timesteps,Grain_number_of_grain_boundaries_for_chunk_timesteps, Grain_number_of_triple_lines_for_chunk_timesteps,Grain_average_misorientation_for_chunk_timesteps,Grain_length_of_triple_lines_for_chunk_timesteps,grain_total_curvature_surface_area_volume_for_chunk_timesteps, GrainBoundary_number_of_triple_line_for_chunk_timesteps,GrainBoundary_misorientation_for_chunk_timesteps,GrainBoundary_length_of_triple_line_for_chunk_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_chunk_timesteps,grain_total_curvature_surface_area_volume_from_face_for_chunk_timesteps,grain_boundary_property_extended_for_chunk_timesteps, EdgeLength_for_chunk_timesteps,InnerAtom_for_chunk_timesteps,GrainBoundary_for_chunk_timesteps,TripleLine_for_chunk_timesteps,QuadraplePoint_for_chunk_timesteps
        
        gc.collect()
