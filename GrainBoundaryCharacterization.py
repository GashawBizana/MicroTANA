import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math
import sys
import random
import time
import concurrent.futures
import multiprocessing as mp
import GrainDataConstructorCurrent
import UnallocatedOrphanAtomGrainIdUpdater


print('usage: GrainBoundaryChracterization.py ./rootFolder/')
#args = sys.argv
rootFolder = 'D:/md_data_Gashaw/'
#args[1]
subDirs = ['AnnealAl161New/Quenched/GradeA_Al161_500']
fname='dump.annealAl161__AtomData_1_165025.cfg'
#file='D:/md_data_Gashaw/AnnealAl161New/Quenched/GradeA_Al161_500/dump.annealAl161__AtomData_1_165025.cfg'#(rootFolder + subDirs + fname)
file='E:/'+fname
Deviation=0 # Grain Width Relaxation
GrainWidth=5 # Average Grain Width
TheresholdDistance= GrainWidth + Deviation # Thereshold Distance for finding nearest neighbours

SeedTestAtom=87606#6596#62841
###################################################
def GrainEdgeLength(FileToAnalize):
    with open(FileToAnalize) as file:
        for line in file:
            if ('Number of particles ' in line):
                NumberOfAtoms_Str=line.split("= ",1)[1]
                NumberOfAtoms=int(NumberOfAtoms_Str[0:-1])
            elif ('H0(1,1)' in line):
                xl_str=line.split("= ",1)[1]
                xl=float(xl_str[0:-1])
            elif ('H0(2,2)' in line):
                yl_str=line.split("= ",1)[1]
                yl=float(yl_str[0:-1])
            elif ('H0(3,3)' in line):
                zl_str=line.split("= ",1)[1]
                zl=float(zl_str[0:-1])
    return(NumberOfAtoms,xl,yl,zl)
################################################
def HeaderFile(file):
    with open(file) as myfile:
        head = [next(myfile) for x in range(24)]
    return(head)
      
    
def ReadGrainData(file):
    GrainData=np.loadtxt(file, skiprows = 24)
    return(GrainData)
        
def EquilidianDistance (X_seed,X,isPeriodic,EdgeLength):

    if isPeriodic:
        dX_NotPeriodic=X-X_seed
        dX_Periodic=dX_NotPeriodic-EdgeLength*((dX_NotPeriodic/EdgeLength).astype(int))
        NeighbourDistance=np.linalg.norm(dX_Periodic,axis=1)
    else:
        dX_NotPeriodic=X-X_seed
        NeighbourDistance=np.linalg.norm(dX_NotPeriodic,axis=1)
    return(NeighbourDistance)
    
def FindNearestNeighbour(GrainId, NeighbourDistance,TheresholdDistance):
    SortedNeighbourDistance=np.sort(NeighbourDistance,kind='mergesort')
    IndexOfSortedNeighbourDistance=np.argsort(NeighbourDistance,kind='mergesort')
    NeigbourList=[]
    i=1
    while (SortedNeighbourDistance[i]<=TheresholdDistance):
        NeigbourList.append(GrainId[IndexOfSortedNeighbourDistance[i],0])
        # if (SortedNeighbourDistance[i]>TheresholdDistance):
        #       break
        i += 1
    # for row,item in enumerate(GrainId):
    #     if (NeighbourDistance[row] <= TheresholdDistance and NeighbourDistance[row]!=0):
    #         NeigbourList.append(item)
    return(NeigbourList)

def FindNumberOfEachNeighbourGrainType(NeighbourList):
    (uniq, freq) =(np.unique(NeighbourList, return_counts=True))
    NumberOfEachNeighbourGrainType=np.column_stack((uniq,freq))[::-1]
    #NumberOfEachNeighbourGrainType=NumberOfEachNeighbourGrainType[::-1]
    return(NumberOfEachNeighbourGrainType)
        
  
def ComputeProportion(NumberOfEachNeighbourGrainType,NeighbourList):
    
    NumberOfNeigbours=len(NeighbourList)
    Proportion=np.array(NumberOfEachNeighbourGrainType[:,[1]])/np.array(NumberOfNeigbours)
    return(Proportion)

def ComputeDiversityIndex(Proportion):
    DiversityIndex=(np.square(Proportion)).sum()
    return(DiversityIndex)

def AtomTopologyClassification(DiversityIndex,NumberOfEachNeighbourGrainType):
    ################### Diversity Index Partition Parameter ###################
    DiversityIndexThrushold_1=1/4
    DiversityIndexThrushold_2=1/12
    DiversityIndexThrushold_3=1/24
    DiversityIndexThrushold_4=1/8
    Partition_GrainToGrainBoundary=0+DiversityIndexThrushold_1 # or 1/2-DiversityIndexThrushold_1
    Partition_GrainBoundaryToTripleLine=1/2+DiversityIndexThrushold_2 # or 2/3-DiversityIndexThrushold_2
    Partition_TripleLineToQuadraplePoint=2/3+DiversityIndexThrushold_3 # or 3/4-DiversityIndexThrushold_3
    Partition_QuadraplePointToChaos=3/4+DiversityIndexThrushold_4 # or 1-DiversityIndexThrushold_4
    ###########################################################################
    if (Partition_GrainToGrainBoundary<=DiversityIndex<=Partition_GrainBoundaryToTripleLine):
        ModifiedGrainId=-1 # Grain Boundary
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1],0]
    elif (Partition_GrainBoundaryToTripleLine<=DiversityIndex<=Partition_TripleLineToQuadraplePoint):
        ModifiedGrainId=-2 # Triple Line
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1,2],0]
    elif (Partition_TripleLineToQuadraplePoint<=DiversityIndex<=Partition_QuadraplePointToChaos):
        ModifiedGrainId=-3 # Quadraple point
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1,2,3],0]
    elif (0<=DiversityIndex<Partition_GrainToGrainBoundary):
        ModifiedGrainId=(NumberOfEachNeighbourGrainType[0,[0]])[0] # Grain
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0],0]
    elif (Partition_QuadraplePointToChaos < DiversityIndex <= 1):
        ModifiedGrainId=-4  # Chaos
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[:,0]

    return(int(ModifiedGrainId),TopologyCreatingNeighbour)
#############################################################

def WholeComputationForOneSeedAtom(Index):
    
    X_seed=X[Index]    
    NeighbourDistance=EquilidianDistance(X_seed,X,isPeriodic,EdgeLength)
    NeighbourList=FindNearestNeighbour(GrainId,NeighbourDistance,TheresholdDistance)
    NumberOfEachNeighbourGrainType=FindNumberOfEachNeighbourGrainType(NeighbourList)
    Proportion=ComputeProportion(NumberOfEachNeighbourGrainType,NeighbourList)

    ModifiedGrainIdOneSeed=int((NumberOfEachNeighbourGrainType[0,[0]])[0])
    TopologyCreatingNeighbourOneSeed=NumberOfEachNeighbourGrainType[[0],0]
    result=[Index,ModifiedGrainIdOneSeed,TopologyCreatingNeighbourOneSeed]
    if (len(Proportion)!=1):
        GINN_SimpsonIndex=1-ComputeDiversityIndex(Proportion)
        ModifiedGrainIdOneSeed,TopologyCreatingNeighbourOneSeed=AtomTopologyClassification(GINN_SimpsonIndex,NumberOfEachNeighbourGrainType)
        result=[Index,ModifiedGrainIdOneSeed,TopologyCreatingNeighbourOneSeed]
    return(result)

ModifiedGrainIdList=[]

def get_ModifiedGrainId(result):
    global ModifiedGrainIdList 
    ModifiedGrainIdList.append(result)
    

            
    
########################################################


NumberOfAtoms,xl,yl,zl=GrainEdgeLength(file)
EdgeLength=np.array([xl,yl,zl])
isPeriodic=True
GrainData=ReadGrainData(file)
X=GrainData[:,[0,1,2]]*EdgeLength
GrainId=GrainData[:,[6]]
UnallocatedOrphanAtomGrainIdUpdater.UpdateGrainIdForUnallocatedOrphanAtoms(GrainId, X,EdgeLength)
if __name__ == '__main__':
    # Index=array=np.arange(0,len(X))
    TopologyCreatingNeighbourList=[]
    startTime = time.time()
    
    pool=mp.Pool(mp.cpu_count())
    
    #X.shape[0]
    for i in range (0,5):
        pool.apply_async(WholeComputationForOneSeedAtom, args=(i,),callback=get_ModifiedGrainId)
    pool.close()
    pool.join()
        # ModifiedGrainIds,TopologyCreatingNeighbours=executor.map(WholeComputationForOneSeedAtom,Index)
        # for ModifiedGrainId in ModifiedGrainIds:
        #     ModifiedGrainIdList.append(ModifiedGrainId)
        #     print(ModifiedGrainId)
        # for TopologyCreatingNeighbour in TopologyCreatingNeighbours:
        #     TopologyCreatingNeighbourList.append(TopologyCreatingNeighbour)
    ModifiedGrainIdListSorted=sorted(ModifiedGrainIdList, key=lambda x: x[0])
    InnerAtom,GrainBoundary,TripleLine,QuadraplePoint=GrainDataConstructorCurrent.OrderedTopologies(ModifiedGrainIdListSorted,X)
    M=[]
    for i in range(0,len(ModifiedGrainIdListSorted)):
        M.append(ModifiedGrainIdListSorted[i][1])
    MM=np.transpose(np.array(M,ndmin=2))
    all_data = np.hstack((GrainData[:,[0,1,2,3,4,5]],GrainId,GrainData[:,[7,8]], MM))
    Header="".join(map(str,HeaderFile(file)))
    f = open('Test2.cfg', 'w')
    f.write(Header)
    f.close()
    with open('Test2.cfg', 'a+') as outfile:
        np.savetxt(outfile,all_data, fmt="%5f")
    executionTime = (time.time() - startTime)
    
    print('Execution time in seconds: ' + str(executionTime))   

##############################Bettter Set For parallel##############################################
# startTime = time.time()
# xl,yl,zl=GrainEdgeLength(file)
# EdgeLength=np.array([xl,yl,zl])
# GrainData=ReadGrainData(file)
# X=GrainData[:,[0,1,2]]*EdgeLength
# GrainId=GrainData[:,[6]]

# ModifiedGrainIdList=[]
# TopologyCreatingNeighbourList=[]

# for X_seed in X:
    
#     ModifiedGrainId,TopologyCreatingNeighbour=WholeComputationForOneSeedAtom(GrainId,X_seed,X)
#     ModifiedGrainIdList.append(ModifiedGrainId)
#     TopologyCreatingNeighbourList.append(TopologyCreatingNeighbour)
    
# executionTime = (time.time() - startTime)
# print('Execution time in seconds: ' + str(executionTime))   



###########################PreviousMethodWithoutParallization#############################################
# xl,yl,zl=GrainEdgeLength(file)
# EdgeLength=np.array([xl,yl,zl])
# GrainData=ReadGrainData(file)
# # x=GrainData[:, [0]]*xl
# # y=GrainData[:, [1]]*yl
# # z=GrainData[:, [2]]*zl
# X=GrainData[:,[0,1,2]]*EdgeLength #[GrainData[:, [0]]*xl,GrainData[:, [1]]*yl, GrainData[:, [2]]*zl]
# GrainId=GrainData[:,[6]]
# startTime = time.time()
# ModifiedGrainIdList=[]
# for SeedAtom in GrainData:
#     X_seed=[SeedAtom[0], SeedAtom[1], SeedAtom[2]]*EdgeLength
#     NeighbourDistance=EquilidianDistance(X_seed,X)
#     NeighbourList=FindNearestNeighbour(GrainId,NeighbourDistance,TheresholdDistance)
#     NumberOfEachNeighbourGrainType=FindNumberOfEachNeighbourGrainType(NeighbourList)
#     Proportion=ComputeProportion(NumberOfEachNeighbourGrainType,NeighbourList)
#     GINN_SimpsonIndex=1-ComputeDiversityIndex(Proportion)
#     ModifiedGrainId,TopologyCreatingNeighbour=AtomTopologyClassification(GINN_SimpsonIndex,NumberOfEachNeighbourGrainType)
#     ModifiedGrainIdList.append(ModifiedGrainId)
# executionTime = (time.time() - startTime)
# print('Execution time in seconds: ' + str(executionTime))

########################################################################################################