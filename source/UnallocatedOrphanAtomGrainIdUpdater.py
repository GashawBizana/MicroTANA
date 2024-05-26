import numpy as np

###### Find the GrainId of the neigbours ####

def FindNumberOfEachNeighbourGrainType(NeighbourList):
    (uniq, freq) =(np.unique(NeighbourList, return_counts=True))
    NumberOfEachNeighbourGrainType=np.column_stack((uniq,freq))[::-1]
    #NumberOfEachNeighbourGrainType=NumberOfEachNeighbourGrainType[::-1]
    return(NumberOfEachNeighbourGrainType) 

def UpdateGrainIdForUnallocatedOrphanAtoms(PeriodicNeighbourList_Index,GrainId,X,EdgeLength,isPeriodic):
    '''
    This function assignes a grain ID for atoms whose grain ID was not assigned by Grade-A 
    The Grain ID for these atoms is assigned to the Grain ID that corrospond to most of its nearby atoms
    '''
    for i in range(0,len(PeriodicNeighbourList_Index)):
        if (GrainId[i]==np.array([-1])):
            PeriodicNeighbourIndex=PeriodicNeighbourList_Index[i]
            NeighbourListGrainId=GrainId[PeriodicNeighbourIndex[:]]
            NumberOfEachGrainType=FindNumberOfEachNeighbourGrainType(NeighbourListGrainId)
            GrainId[i]=(NumberOfEachGrainType[0,[0]])[0]
            
        

