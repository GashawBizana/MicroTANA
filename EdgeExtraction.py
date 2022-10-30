import numpy as np
import open3d as o3d
xyz=[]
for i in range (0, len(X)):

    if GrainId[i] ==[0]:
        
        xyz.append(X[i])
xyz=GrainBoundary1[0]
        
pcd= o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(xyz)
# o3d.visualization.draw_geometries([pcd])
# kdt=pcd.o3d.geometry.KDTreeSearchParamKNN(knn= 30)

Cov=np.asarray(o3d.geometry.PointCloud.estimate_point_covariances(pcd,
    search_param = o3d.geometry.KDTreeSearchParamKNN(knn=10)))

eig=np.linalg.eig(Cov[:])[0]

Var=np.min(eig,axis=1)/np.sum(eig[:],axis=1)
T=0.1
Index=np.where(Var>0.2)[0]
Test2=np.take(xyz,Index,axis=0)
pcd= o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(Test2)
o3d.visualization.draw_geometries([pcd])