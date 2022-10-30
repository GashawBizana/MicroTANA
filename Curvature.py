import numpy as np

from NormalEstimation import GaussianWeightFunction
## Gausssian weigh function centered around a sample point x and a neighbour q
# def GaussianWeightFunction(x,Q,h):
#     Theta_Q=[]
#     for i in range(0,len(Q)):
#         Theta_q= np.exp(-np.square(np.linalg.norm(x-Q[i]))/(h*h))
#         Theta_Q.append(Theta_q)
#     return(Theta_Q)


def VectorField(Normal_Q,Theta_Q):
    Inner=[]
    for i in range(0,len(Theta_Q)):
        Normalq_iTimesThetaq_i=Normal_Q[i]*Theta_Q[i]
        Inner.append(Normalq_iTimesThetaq_i)
    N_x=np.sum(Inner,axis = 0)
    n_x=N_x/np.linalg.norm(N_x)
    
    return(n_x,N_x,Inner)
    
def g_of_x(x,Theta_Q,Q,n_x):
    Inner=[]
    for i in range(0,len(Q)):
        InnerOfg_of_x_i=(2*Theta_Q[i])*(np.dot(x-Q[i],n_x))
        Inner.append(InnerOfg_of_x_i)
    g_x=np.sum(Inner,axis = 0)
    
    return(g_x)
        
def XForCurvature(NormalVector,PeriodicNeighbourList_Index,X,h,EdgeLength,Periodic):
    X_forCurvature=[]
    for j in range(0,len(PeriodicNeighbourList_Index)):
        PeriodicNeighbourIndex=PeriodicNeighbourList_Index[j]
        #PeriodicNeighbourIndex.remove(j)
        Q=X[PeriodicNeighbourIndex[:]]
        x=X[j]
        error=1
        while (error > 0.00001):
            xtemp=x
            Normal_Q=NormalVector[PeriodicNeighbourIndex[:]]
            Theta_Q=GaussianWeightFunction(x,Q,h,EdgeLength,Periodic)
            n_x,N_x,Inner=VectorField(Normal_Q,Theta_Q)
            g_x= g_of_x(x,Theta_Q,Q,n_x)
            x=x+(g_x/np.sum(Theta_Q,axis = 0))*n_x
            error=np.linalg.norm(x-xtemp)
        X_forCurvature.append(x)
            
# x=X1[0]
# PeriodicNeighbourIndex=PeriodicNeighbourList_Index[0]
# Q=X1[PeriodicNeighbourIndex[:]]
# Theta_Q=GaussianWeightFunction(x,Q,h,EdgeLength,Periodic)
# n_x,N_x,Inner=VectorField(Normal_Q,Theta_Q)

# X_forCurvature=[]
# x=XX[0]
# error=1
# while (error > 0.00001):
#     xtemp=x
#     PeriodicNeighbourIndex=PeriodicNeighbourList_Index[0]
#     Q=XX[PeriodicNeighbourIndex[:]]
#     Normal_Q=NormalVector[PeriodicNeighbourIndex[:]]
#     Theta_Q=GaussianWeightFunction(x,Q,h,EdgeLength,Periodic)
#     n_x,N_x,Inner=VectorField(Normal_Q,Theta_Q)
#     print(n_x)
#     g_x= g_of_x(x,Theta_Q,Q,n_x)
#     print(g_x)
#     x=x+(g_x/np.sum(Theta_Q,axis = 0))*n_x
#     error=np.linalg.norm(x-xtemp)
#     print(error)
#     X_forCurvature.append(x)