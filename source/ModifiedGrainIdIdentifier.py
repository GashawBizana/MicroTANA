
import numpy as np



######## This file contain F=functions used to determine Modified GrainId of one seed atom given its neighbourList ###

###### Find the GrainId of the neigbours ####

def FindNumberOfEachNeighbourGrainType(NeighbourList):
    (uniq, freq) =(np.unique(NeighbourList, return_counts=True))
    NumberOfEachNeighbourGrainType=np.column_stack((uniq,freq))[::-1]
    #NumberOfEachNeighbourGrainType=NumberOfEachNeighbourGrainType[::-1]
    return(NumberOfEachNeighbourGrainType)


def ComputeProportion(NumberOfEachNeighbourGrainType,NeighbourList):
     ######## This function compute the proportion of each neighbour GrainId in the neighbour list ###    
    NumberOfNeigbours=len(NeighbourList)
    Proportion=np.array(NumberOfEachNeighbourGrainType[:,[1]])/np.array(NumberOfNeigbours)
    return(Proportion)


def ComputeDiversityIndex(Proportion,NumberOfEachGrainType):

    ############ This function computes the Simpson Diversity Index  ###

    ### Here we are using the simpson diversity index
    numberofspecies0=ComputeDiversityIndexSimpson(Proportion)
    eveness0=ComputeDiversityIndexShannon(Proportion)/np.log(numberofspecies0)
    
    
    DiversityIndex=numberofspecies0
    
     # increase 0.6 to consider more of the higher order topological defects (TL and QP) as lower order topological defect
    if 0.9 < eveness0 and (numberofspecies0==3 or numberofspecies0==4) :
        if numberofspecies0==3:
        
            numberofspecies1=numberofspecies0-1
            DiversityIndex=numberofspecies1
            
        
        elif numberofspecies0==4:
            
            numberofspecies1=numberofspecies0-2
            DiversityIndex=numberofspecies1
           
    
    return(DiversityIndex)  
        
        
        
    
    
    

def ComputeDiversityIndexSimpson(Proportion):
    # An alternative way to compute diversity index and its inverse from 
    DiversityIndex=(np.square(Proportion)).sum()
    DiversityIndex=1/DiversityIndex
    
    DiversityIndex=len(Proportion)
    return(DiversityIndex)

def ComputeDiversityIndexShannon(Proportion):
    DiversityIndex=-(Proportion*(np.log(Proportion))).sum()
    
    return(DiversityIndex)


########## Classfies which Microstructure network entity an atom belongs to based on Diversity index ###
def AtomTopologyClassification(DiversityIndex,NumberOfEachNeighbourGrainType):
    ################### Diversity Index Partition Parameter ###################
    # DiversityIndexThrushold_1=0 #1/4
    # DiversityIndexThrushold_2=1/12
    # DiversityIndexThrushold_3=1/24
    # DiversityIndexThrushold_4=1/8
    # Partition_GrainToGrainBoundary=0+DiversityIndexThrushold_1 # or 1/2-DiversityIndexThrushold_1
    # Partition_GrainBoundaryToTripleLine=1/2+DiversityIndexThrushold_2 # or 2/3-DiversityIndexThrushold_2
    # Partition_TripleLineToQuadraplePoint=2/3+DiversityIndexThrushold_3 # or 3/4-DiversityIndexThrushold_3
    # Partition_QuadraplePointToChaos=3/4+DiversityIndexThrushold_4 # or 1-DiversityIndexThrushold_4
    

    Partition_GrainToGrainBoundary=1.0 # or 1/2-DiversityIndexThrushold_1
    Partition_GrainBoundaryToTripleLine=2.01 # or 2/3-DiversityIndexThrushold_2
    Partition_TripleLineToQuadraplePoint=3.01 # or 3/4-DiversityIndexThrushold_3
    Partition_QuadraplePointToChaos=7 # or 1-DiversityIndexThrushold_4
    
    ###########################################################################
    if (Partition_GrainToGrainBoundary<DiversityIndex<Partition_GrainBoundaryToTripleLine):
        ModifiedGrainId=-1 # Grain Boundary
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1],0]
    elif (Partition_GrainBoundaryToTripleLine<=DiversityIndex<Partition_TripleLineToQuadraplePoint):
        ModifiedGrainId=-2 # Triple Line
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1,2],0]
    elif (Partition_TripleLineToQuadraplePoint<=DiversityIndex<=Partition_QuadraplePointToChaos):
        ModifiedGrainId=-3 # Quadraple point
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1,2,3],0]
    elif (1<=DiversityIndex<Partition_GrainToGrainBoundary):
        ModifiedGrainId (NumberOfEachNeighbourGrainType[0,[0]])[0] # Grain
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0],0]
    elif (Partition_QuadraplePointToChaos < DiversityIndex):
        ModifiedGrainId=-4  # Chaos
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[:,0]
    # else:
    #     ModifiedGrainId=(NumberOfEachNeighbourGrainType[0,[0]])[0] # Grain
    #     TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0],0]
    return(int(ModifiedGrainId),TopologyCreatingNeighbour)





########### Identifies the Modified GrainId for all atoms

def ModifiedGrainIdComputation(PeriodicNeighbourList_Index,GrainId,X,EdgeLength,isPeriodic):
    '''
    Given a list of neighbours for each atom and their corrosponding grain ID, this function assignes a 
    "Modified Grain ID" which inicates wheather an atoms is part of an inside of a grain, a grain boundary, 
    a triple junction or a quadraple junction.
    '''
    ModifiedGrainIdList=[]
    for i in range(0,len(PeriodicNeighbourList_Index)):
        PeriodicNeighbourIndex=PeriodicNeighbourList_Index[i]
        NeighbourListGrainId=GrainId[PeriodicNeighbourIndex[:]]
        NumberOfEachGrainType=FindNumberOfEachNeighbourGrainType(NeighbourListGrainId)
        ProportionOfEachGrainType=ComputeProportion(NumberOfEachGrainType,NeighbourListGrainId)
        ModifiedGrainIdOfOneSeed=int((NumberOfEachGrainType[0,[0]])[0])
        TopologyCreatingNeighboursOfOneSeed=NumberOfEachGrainType[[0],0]
        result=[i,ModifiedGrainIdOfOneSeed,TopologyCreatingNeighboursOfOneSeed]
        if (len(ProportionOfEachGrainType)!=1):
            #GINN_SimpsonIndex=ComputeDiversityIndexSimpson(ProportionOfEachGrainType)
            GINN_SimpsonIndex=ComputeDiversityIndex(ProportionOfEachGrainType,NumberOfEachGrainType)
            ModifiedGrainIdOfOneSeed,TopologyCreatingNeighboursOfOneSeed=AtomTopologyClassification(GINN_SimpsonIndex,NumberOfEachGrainType)
        '''
        This was an initial  attempt to compute per atom curvature using the method by Mason 
        by calling the function at the bottom as:
        CurvatureMason(X_seed,ModifiedGrainIdOfOneSeed,GrainId_seed,NeighbourListPosition,NeighbourListGrainId,EdgeLength,isPeriodic)
        '''
        result=[i,ModifiedGrainIdOfOneSeed,TopologyCreatingNeighboursOfOneSeed]
            
        
        ModifiedGrainIdList.append(result)
    return(ModifiedGrainIdList)





















