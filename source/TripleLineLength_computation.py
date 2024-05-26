import numpy as np
from sklearn.cluster import DBSCAN
from scipy.spatial import ConvexHull
from scipy.spatial.distance import cdist

import scipy.spatial as spatial

from scipy.interpolate import splprep, splev
import copy


# def UnWrapCluster(pcd,EdgeLength):
#     # NumberOfPoints=len(pcd.points)
#     labels = np.array(pcd.cluster_dbscan(eps=5, min_points=1))
#     max_label = labels.max()
#     id_big=np.where(labels == 0)[0]
#     pcd_big = pcd.select_by_index(id_big)
#     pcd_all=pcd_big
        
    
#     for i in range(1,max_label+1):
        
#         id_i=np.where(labels == i)[0]
#         pcd_i = pcd.select_by_index(id_i)
#         diff_center=np.rint((pcd_big.get_center()-pcd_i.get_center())/EdgeLength)
#         pcd_i_to_big=pcd_i.translate(diff_center*EdgeLength,relative=True)
#         pcd_all += pcd_i_to_big
        
#     # labels_2 = np.array(pcd_all.cluster_dbscan(eps=5, min_points=1))    
#     # id_big_2=np.where(labels_2 == 0)[0]
#     # pcd_big_2 = pcd_all.select_by_index(id_big_2)
    
#     # for k in range(0,labels_2.max()+1):
#     #     id_2_k=np.where(labels_2  == k)[0]
#     #     pcd_2_k = pcd_all.select_by_index(id_2_k)
        
#     #     if len(pcd_2_k.points ) > len(pcd_big_2.points):
#     #         pcd_big_2=pcd_all.select_by_index(id_2_k)
        
        
#     # pcd_as_array=np.asarray(pcd_big_2.points)
    
#     pcd_as_array=np.asarray(pcd_all.points)
#     return(pcd_as_array)

def UnWrapCluster(pcd,EdgeLength):

    model = DBSCAN(eps=5, min_samples=1)
    model.fit_predict(pcd)
    labels=model.labels_ 
    max_label=labels.max()
    id_big=np.where(labels == max_label)[0]
    pcd_big = pcd[id_big]
    
    d=[]
    for i in range(0,max_label+1):
        id_i=np.where(labels == i)[0]
        if id_i.size==0:
            continue
        pcd_i = pcd[id_i]
        diff_center=np.rint((np.mean(pcd_big,axis=0)-np.mean(pcd_i,axis=0))/EdgeLength)
        d.append(diff_center)
    pu_ac=copy.deepcopy(pcd)
    
    pu_ac=pu_ac+(np.array(d)[labels])*EdgeLength
    return(pu_ac)


def TripleLineStartOrEnd(thinned_points):
    
    points=thinned_points
    hull = points#ConvexHull(points)

    # Extract the points forming the hull
    hullpoints = hull #points[hull.vertices,:]
    
    # Naive way of finding the best pair in O(H^2) time if H is number of points on
    # hull
    hdist = cdist(hullpoints, hullpoints, metric='euclidean')
    
    # Get the farthest apart points
    bestpair = np.unravel_index(hdist.argmax(), hdist.shape)
    
    P = np.array([hullpoints[bestpair[0]],hullpoints[bestpair[1]]])
    
    # Now we have a problem
    
    # uncomment below if needed
    # while len(P)<1:
    #     # was one
        
     
    #   distance_to_P        = cdist(points, P)
    #   minimum_to_each_of_P = np.min(distance_to_P, axis=1)
    #   best_new_point_idx   = np.argmax(minimum_to_each_of_P)
    #   best_new_point = np.expand_dims(points[best_new_point_idx,:],0)
    #   P = np.append(P,best_new_point,axis=0)
    
    return (P[0],np.max(hdist))



def thin_line(points, point_cloud_thickness=6, iterations=1,sample_points=0):
    
    if sample_points != 0:
        points = points[:sample_points]
    
    # Sort points into KDTree for nearest neighbors computation later
    point_tree = spatial.cKDTree(points)

    # Empty array for transformed points
    new_points = []
    # Empty array for regression lines corresponding ^^ points
    regression_lines = []
   
    for point in point_tree.data:
        # Get list of points within specified radius {point_cloud_thickness}
        
        points_in_radius = point_tree.data[point_tree.query_ball_point(point, point_cloud_thickness)]
        

        # Get mean of points within radius
       
        data_mean = points_in_radius.mean(axis=0)

        # Calulate 3D regression line/principal component in point form with 2 coordinates
        uu, dd, vv = np.linalg.svd((points_in_radius - data_mean)**1)
        linepts = vv[0] * np.mgrid[-1:1:2j][:, np.newaxis]
        linepts += data_mean
        regression_lines.append(list(linepts))
        

        # Project original point onto 3D regression line
       
        ap = point - linepts[0]
        ab = linepts[1] - linepts[0]
        point_moved = linepts[0] + np.dot(ap,ab) / np.dot(ab,ab) * ab

        new_points.append(list(point_moved))
    return np.array(new_points), regression_lines



def distance(P1, P2):
    """
    This function computes the distance between 2 points defined by
     P1 = (x1,y1) and P2 = (x2,y2) 
    """

    return ((P1[0] - P2[0])**2 + (P1[1] - P2[1])**2 + (P1[2] - P2[2])**2) ** 0.5


def optimized_path(coords, start=None):
    """
    This function finds the nearest point to a point
    coords should be a list in this format coords = [ [x1, y1], [x2, y2] , ...] 
    
    """
    if start is None:
        start = coords[0]
    pass_by = coords
    path = [start]
    pass_by.remove(start)
    while pass_by:
        nearest = min(pass_by, key=lambda x: distance(path[-1], x))
        path.append(nearest)
        pass_by.remove(nearest)
    return np.array(path)

def curve_length(curve):
    """ sum of Euclidean distances between points """
    return np.sum(np.sqrt(np.sum((curve[:-1] - curve[1:])**2,axis=0)))

def curve_length1(curve):
    
    s=0
    for i in range(0, len(curve)-1):
        
        if np.linalg.norm(curve[i]-curve[i+1]) < 5:
            #7 
        
            s=s+np.linalg.norm(curve[i]-curve[i+1])
    return(s)
    

def interpolate_otimized_path(triple_line_path):
    # very small error added to avoid repeated values
    n1=np.random.normal(0, 0.00000001, len(triple_line_path))
    n2=np.random.normal(0, 0.00000001, len(triple_line_path))
    n3=np.random.normal(0, 0.00000001, len(triple_line_path))
    
    noise=np.array( [n1,n2,n3])
    
    triple_line_path=triple_line_path.T+noise
    if len(triple_line_path) <=3:
        new_points = triple_line_path
    elif len(triple_line_path) >3:
        
        tck, u = splprep(triple_line_path, s=4)
        #s=1
        new_points = splev(u, tck)
    ret=np.array(new_points)
    return(ret.T)
    

def TripleLineLength(triple_line_points,EdgeLength):
    
    triple_line_points_l=UnWrapCluster(triple_line_points, EdgeLength)
    thinned_points, regression_lines=thin_line(triple_line_points_l, point_cloud_thickness=7)
    # 4
    P_start,lam_i=TripleLineStartOrEnd(thinned_points)
    
    triple_line_path=optimized_path(thinned_points.tolist(),P_start.tolist())
    
    ### comment for only one iteration
    # P_start=TripleLineStartOrEnd(triple_line_path)
    
    # triple_line_path=optimized_path(triple_line_path.tolist(),P_start.tolist())
    ##
    interpolated_triple_line_path=interpolate_otimized_path(triple_line_path)
    
    triple_line_length= curve_length1(interpolated_triple_line_path)
    
    return(triple_line_length,interpolated_triple_line_path,lam_i)


def TripleLineLength_For_all(TripleLine,EdgeLength):
    TripleLineLength_for_all=[]
    for i in range(0,len(TripleLine)):
        
        triple_line_ID_i= TripleLine[i][0]
        
        triple_line_points_i=np.array(TripleLine[i][1:])
        
        triple_line_length_i,interpolated_triple_line_path_i,lam_i=TripleLineLength(triple_line_points_i,EdgeLength)
         
        TripleLineLength_for_all.append([triple_line_ID_i,triple_line_length_i,interpolated_triple_line_path_i,lam_i])    
    return(TripleLineLength_for_all)
    
    
# fig2 = plt.figure(2)
# ax3d = fig2.add_subplot(111, projection='3d')
# for i in range(0,len(TripleLineLength_for_all)):
#     if 0 in TripleLineLength_for_all[i][0]:
#       TL=TripleLineLength_for_all[i][2]
  
#       ax3d.plot(TL[:,0], TL[:,1], TL[:,2], 'bo-')
# plt.show()


# fig3 = plt.figure(3)
# ax3d = fig3.add_subplot(111, projection='3d')
# for i in range(0,len(TripleLine)):
#     if 0 in TripleLine[i][0]:
#       TL=np.array(TripleLine[i][1:])
  
#       ax3d.plot(TL[:,0], TL[:,1], TL[:,2], 'bo')
# plt.show()

# TripleLine=triple_line_length_for_all_timesteps[0]
# fig3 = plt.figure(3)
# ax3d = fig3.add_subplot(111, projection='3d')
# TL=TripleLine[0][2]
  
# ax3d.plot(TL[:,0], TL[:,1], TL[:,2], 'bo-')
# plt.show()
'''
j=0
tlll=triple_line_length_for_all_timesteps[0]
gm=grain_property_for_all_timesteps[0][j][1]
fig2 = plt.figure(2)
ax3d = fig2.add_subplot(111, projection='3d')
for i in range(0,len(tlll)):
    if j in tlll[i][0]:
      TL=tlll[i][2]
      
      ax3d.plot(TL[:,0], TL[:,1], TL[:,2],"o--", markersize=5)
plt.show() 
# gm.show(viewer='gl')
'''