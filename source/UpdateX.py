import numpy as np
from Helper_functions_for_main import PeriodicDistance, PeriodicCenterOfMass

def update_x(X,ModifiedGrainIdListSorted,GrainId,EdgeLength,isPeriodic,PeriodicNeighbourList_Index):
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
    return np.array(Xnew)