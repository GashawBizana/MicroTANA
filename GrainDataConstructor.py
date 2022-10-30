
######## This Module Creats a suitable data structure of Microstricture Network entities

###  Groups all atoms that belong to a cetain GB, TL, QP or IA together without duplication ###
def TopologyGrouper (Topology):
    new_Topology = []
    for ParentGrains, AtomPositionOfTopology in Topology:
        for index, item in enumerate(new_Topology):
            if sorted(item[0]) == sorted(ParentGrains):
                # if AtomPositionOfTopology not in item:
                new_Topology[index].append(AtomPositionOfTopology)
                break
        else:
            new_Topology.append([ParentGrains, AtomPositionOfTopology])
    new_Topology=sorted(new_Topology, key=lambda l: int(l[0][0]),reverse=False)
    return(new_Topology)

# def TopologyGrouper1 (Topology):
#     new_Topology = []
#     for i in range(0,len(Topology)):
#         for index, item in enumerate(new_Topology):
#             if sorted(new_Topology[j]) == sorted(Topology[i][0]):
#                 if Topology[i][1] not in new_
    

##### Computs the numbe of Grains that are in our microstructure ###
def NumberOfGrainCalculator(ModifiedGrainIdListSorted):
    NumberOfGrains=int(max(ModifiedGrainIdListSorted, key=lambda x: x[1])[1])+1
    return(NumberOfGrains)

########## Creats Structured data for each each kind of network entity ######

def GrainInfoConstructor(GrainList,X):
    InnerAtom=[]
    GrainBoundary=[]
    TripleLine=[]
    QuadraplePoint=[]
    Curvature=[]
    for i in range(0,len(GrainList)-1):
        ModifiedGrainId=(GrainList[i])[1]
        TopologyCreatingList=(GrainList[i])[2]
        AtomPosition=X[int((GrainList[i])[0])].tolist()
        CurvatureOfEachAtom=(GrainList[i])[3]
        
        if (ModifiedGrainId>=0):
            
            
             InnerAtom.append([[int(TopologyCreatingList[0])],AtomPosition])
            
        elif (ModifiedGrainId ==-1):
            
            GrainBoundary.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[1])],AtomPosition])
            # GrainBoundary.append([[int(TopologyCreatingList[1]),int(TopologyCreatingList[0])],AtomPosition])
            
            Curvature.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[1])],[CurvatureOfEachAtom]])
            
            # InnerAtom.append([[int(TopologyCreatingList[0])],AtomPosition])
            # InnerAtom.append([[int(TopologyCreatingList[1])],AtomPosition])

        elif (ModifiedGrainId ==-2):
            TripleLine.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[1]),int(TopologyCreatingList[2])],AtomPosition])
            
            
            # GrainBoundary.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[1])],AtomPosition])
            # GrainBoundary.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[2])],AtomPosition])
            # GrainBoundary.append([[int(TopologyCreatingList[1]),int(TopologyCreatingList[2])],AtomPosition])
            
            # InnerAtom.append([[int(TopologyCreatingList[0])],AtomPosition])
            # InnerAtom.append([[int(TopologyCreatingList[1])],AtomPosition])
            # InnerAtom.append([[int(TopologyCreatingList[2])],AtomPosition])
            
        elif (ModifiedGrainId ==-3):
            QuadraplePoint.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[1]),int(TopologyCreatingList[2]),int(TopologyCreatingList[3])],AtomPosition])
            
            TripleLine.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[1]),int(TopologyCreatingList[2])],AtomPosition])
            TripleLine.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[1]),int(TopologyCreatingList[3])],AtomPosition])
            TripleLine.append([[int(TopologyCreatingList[1]),int(TopologyCreatingList[2]),int(TopologyCreatingList[3])],AtomPosition])
            TripleLine.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[2]),int(TopologyCreatingList[3])],AtomPosition])
            
            # GrainBoundary.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[1])],AtomPosition])
            # GrainBoundary.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[2])],AtomPosition])
            # GrainBoundary.append([[int(TopologyCreatingList[0]),int(TopologyCreatingList[3])],AtomPosition])
            # GrainBoundary.append([[int(TopologyCreatingList[1]),int(TopologyCreatingList[2])],AtomPosition])
            # GrainBoundary.append([[int(TopologyCreatingList[1]),int(TopologyCreatingList[3])],AtomPosition])
            # GrainBoundary.append([[int(TopologyCreatingList[2]),int(TopologyCreatingList[3])],AtomPosition])
            
            # InnerAtom.append([[int(TopologyCreatingList[0])],AtomPosition])
            # InnerAtom.append([[int(TopologyCreatingList[1])],AtomPosition])
            # InnerAtom.append([[int(TopologyCreatingList[2])],AtomPosition])
            # InnerAtom.append([[int(TopologyCreatingList[3])],AtomPosition])
            
    return(InnerAtom,GrainBoundary,TripleLine,QuadraplePoint,Curvature)

#######   Orders the the structured network entity data according to the Id of the first grain ##
def OrderedTopologies(ModifiedGrainIdListSorted,X):
    #NumberOfGrains=NumberOfGrainCalculator(ModifiedGrainIdListSorted)
    InnerAtom,GrainBoundary,TripleLine,QuadraplePoint,Curvature=GrainInfoConstructor(ModifiedGrainIdListSorted,X)
    new_InnerAtom=TopologyGrouper(InnerAtom)
    new_GrainBoundary=TopologyGrouper(GrainBoundary)
    new_Curvature=TopologyGrouper(Curvature)
    new_TripleLine=TopologyGrouper(TripleLine)
    new_QuadraplePoint=TopologyGrouper(QuadraplePoint)
    return(new_InnerAtom,new_GrainBoundary,new_TripleLine,new_QuadraplePoint,new_Curvature)



        

