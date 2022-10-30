
import numpy as np
import math


######## Functions used to determine Modified GrainId of one seed atom given its neighbourList ###

###### Find the GrainId of the neigbours ####

def FindNumberOfEachNeighbourGrainType(NeighbourList):
    (uniq, freq) =(np.unique(NeighbourList, return_counts=True))
    NumberOfEachNeighbourGrainType=np.column_stack((uniq,freq))[::-1]
    #NumberOfEachNeighbourGrainType=NumberOfEachNeighbourGrainType[::-1]
    return(NumberOfEachNeighbourGrainType)

######## Compute the proportion of each neighbour GrainId in the neighbour list ###
def ComputeProportion(NumberOfEachNeighbourGrainType,NeighbourList):
    
    NumberOfNeigbours=len(NeighbourList)
    Proportion=np.array(NumberOfEachNeighbourGrainType[:,[1]])/np.array(NumberOfNeigbours)
    return(Proportion)

############  Computes the Simpson Diversity Index  ###

### Here we are using the simpson diversity index
def ComputeDiversityIndex(Proportion,NumberOfEachGrainType):
    numberofspecies0=len(Proportion)
    eveness0=ComputeDiversityIndexShannon(Proportion)/np.log(numberofspecies0)
    
    
    DiversityIndex=numberofspecies0
    NumberOfEachGrainType_m=NumberOfEachGrainType
        
    if 0.7 < eveness0 and (numberofspecies0==3 or numberofspecies0==4) :
        if numberofspecies0==3:
        
            numberofspecies1=numberofspecies0-1
            DiversityIndex=numberofspecies1
            NumberOfEachGrainType_m=NumberOfEachGrainType[:numberofspecies1]
        
        elif numberofspecies0==4:
            
            numberofspecies1=numberofspecies0-2
            DiversityIndex=numberofspecies1
            NumberOfEachGrainType_m=NumberOfEachGrainType[:numberofspecies1]
    
    return(DiversityIndex)  
        
        
        
    
    
    

def ComputeDiversityIndexSimpson(Proportion):
    DiversityIndex=(np.square(Proportion)).sum()
    DiversityIndex=1/DiversityIndex
    
    DiversityIndex=len(Proportion)
    return(DiversityIndex)

def ComputeDiversityIndexShannon(Proportion):
    DiversityIndex=-(Proportion*(np.log(Proportion))).sum()
    
    return(DiversityIndex)


########## Classfies which Microstructure network entity an atom belongs to based on DI ###
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
    if DiversityIndex==2:#(Partition_GrainToGrainBoundary<DiversityIndex<Partition_GrainBoundaryToTripleLine):
        ModifiedGrainId=-1 # Grain Boundary
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1],0]
    elif DiversityIndex==3:# (Partition_GrainBoundaryToTripleLine<=DiversityIndex<Partition_TripleLineToQuadraplePoint):
        ModifiedGrainId=-2 # Triple Line
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1,2],0]
    elif DiversityIndex>=4:#(Partition_TripleLineToQuadraplePoint<=DiversityIndex<=Partition_QuadraplePointToChaos):
        ModifiedGrainId=-3 # Quadraple point
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0,1,2,3],0]
    elif DiversityIndex==1:#(1<=DiversityIndex<Partition_GrainToGrainBoundary):
        ModifiedGrainId=(NumberOfEachNeighbourGrainType[0,[0]])[0] # Grain
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0],0]
    elif DiversityIndex>=7:#(Partition_QuadraplePointToChaos < DiversityIndex):
        ModifiedGrainId=-4  # Chaos
        TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[:,0]
    # else:
    #     ModifiedGrainId=(NumberOfEachNeighbourGrainType[0,[0]])[0] # Grain
    #     TopologyCreatingNeighbour=NumberOfEachNeighbourGrainType[[0],0]
    return(int(ModifiedGrainId),TopologyCreatingNeighbour)





########### Identifies the Modified GrainId for all atoms

def ModifiedGrainIdComputation(PeriodicNeighbourList_Index,GrainId,X,EdgeLength,isPeriodic):

    ModifiedGrainIdList=[]
    for i in range(0,len(PeriodicNeighbourList_Index)):
        PeriodicNeighbourIndex=PeriodicNeighbourList_Index[i]
        NeighbourListGrainId=GrainId[PeriodicNeighbourIndex[:]]
        NeighbourListPosition=X[PeriodicNeighbourIndex[:]]
        X_seed=X[i]
        GrainId_seed=GrainId[i]
        NumberOfEachGrainType=FindNumberOfEachNeighbourGrainType(NeighbourListGrainId)
        ProportionOfEachGrainType=ComputeProportion(NumberOfEachGrainType,NeighbourListGrainId)
        ModifiedGrainIdOfOneSeed=int((NumberOfEachGrainType[0,[0]])[0])
        TopologyCreatingNeighboursOfOneSeed=NumberOfEachGrainType[[0],0]
        result=[i,ModifiedGrainIdOfOneSeed,TopologyCreatingNeighboursOfOneSeed]
        if (len(ProportionOfEachGrainType)!=1):
            GINN_SimpsonIndex=ComputeDiversityIndexSimpson(ProportionOfEachGrainType)
            # GINN_SimpsonIndex=ComputeDiversityIndex(ProportionOfEachGrainType,NumberOfEachGrainType)
            ModifiedGrainIdOfOneSeed,TopologyCreatingNeighboursOfOneSeed=AtomTopologyClassification(GINN_SimpsonIndex,NumberOfEachGrainType)
        MeanCurvature=0# CurvatureMason(X_seed,ModifiedGrainIdOfOneSeed,GrainId_seed,NeighbourListPosition,NeighbourListGrainId,EdgeLength,isPeriodic)
        result=[i,ModifiedGrainIdOfOneSeed,TopologyCreatingNeighboursOfOneSeed,MeanCurvature]
            
        
        ModifiedGrainIdList.append(result)
    return(ModifiedGrainIdList)

# Computes the equildian distance between two points -- Same funiction as one in UnallocatedOrphanAtomGrianIdUpdater
def PeriodicDistance(X1,X2,EdgeLength,isPeriodic):
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
    
    d=np.sqrt(np.square(dx)+np.square(dy)+np.square(dz))
    return(d)


# Computs the Mean Curvature of an atom in the GrainBoundary using Mason's paper

def CurvatureMason(x_seed,ModifiedGrainId_seed,GrainId_seed,NeighbourListPosition,NeighbourListGrainId,EdgeLength,isPeriodic):
    if (ModifiedGrainId_seed == -1):
        
        NearestNeighbourListPositionExcludingX_seed=[]
        NeighbourListGrainIdExcludingX_seed=[]
        for j in range(0,len(NeighbourListPosition)):
            if ((x_seed != NeighbourListPosition[j]).any()):
               NearestNeighbourListPositionExcludingX_seed.append(NeighbourListPosition[j]) 
               NeighbourListGrainIdExcludingX_seed.append(NeighbourListGrainId[j])
        Sigma=np.mean(np.std(NearestNeighbourListPositionExcludingX_seed, axis=0))
        D_array=[]
        for i in range(0,len(NeighbourListGrainIdExcludingX_seed)):
            d_ij=PeriodicDistance(x_seed,NearestNeighbourListPositionExcludingX_seed[i],EdgeLength,isPeriodic)
            F_ij=np.exp((np.square(d_ij))/(2*Sigma*Sigma))
            if (F_ij>=0.001): # 0.001 is C in the pape which is usually much less than 1 and we get to choose it
                if (GrainId_seed==NeighbourListGrainIdExcludingX_seed[i]): 
                    Kron_ij=1
                elif(GrainId_seed!=NeighbourListGrainIdExcludingX_seed[i]): 
                    Kron_ij=0
                D_ij=F_ij*(2*Kron_ij-1)
            elif(F_ij<0.001):
                D_ij=0
            D_array.append(D_ij)
            D=np.sum(D_array)
            ####### These values depend on the atomic system and are psuedo elements to the voxel properties in
            # Mason's paper 
            a=4.046  # lattice constant
            beta=0.5 #  A constant implying grainboundary exist between two adjecent grains
            lamda=0.707*a # spacing between atoms or Voxel  spacing  in Mason
            alpha=a**(3/4) # Volume occupied by an atom or Voxel volume 
            H=((np.exp(  ((beta**2)*(lamda**2) ) / (2*(Sigma**2)))) *            
               ( ((alpha*D)/(4*math.pi*(Sigma**2)) )  - ( ( (math.sqrt(math.pi))/(math.sqrt(2)* Sigma) ) * math.erf( (beta*lamda)/(math.sqrt(2*Sigma)) ) ) )) 
    
    else:
        H=0
    return(H)