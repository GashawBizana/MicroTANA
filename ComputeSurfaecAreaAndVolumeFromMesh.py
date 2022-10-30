import numpy as np
import open3d as o3d
import trimesh
import copy
xyz=[]
label0=[]
label1=[]
for i in range (0, len(X)):

    if GrainId[i] ==[0]:
        
        xyz.append(X[i])
        

label1=np.array(label1)

radii=[2]

# ii=
# xyz=np.array(TripleLine[ii][1:])

pcdWrapped = o3d.geometry.PointCloud()
pcdWrapped.points = o3d.utility.Vector3dVector(xyz)

o3d.visualization.draw_geometries([pcdWrapped])

def UnWrapCluster(pcd,EdgeLength):
    NumberOfPoints=len(pcd.points)
    labels = np.array(pcd.cluster_dbscan(eps=5, min_points=1))
    max_label = labels.max()
    id_big=np.where(labels == 0)[0]
    pcd_big = pcd.select_by_index(id_big)
    
    for i in range(0,max_label+1):
        id_i=np.where(labels == i)[0]
        pcd_i = pcd.select_by_index(id_i)
        diff_center=np.rint((pcd_big.get_center()-pcd_i.get_center())/EdgeLength)
        pcd_i_to_big=pcd_i.translate(diff_center*EdgeLength,relative=True)
        pcd_big += pcd_i_to_big
    return(pcd_big)
        
        
        
        
        
        
    
    

# cl, ind = pcdWrapped.remove_radius_outlier(nb_points=int(len(xyz)/3),radius=EdgeLength[0]/2)
# pcd = pcd.select_by_index(ind)

pcd=UnWrapCluster(pcdWrapped,EdgeLength)

pcd.estimate_normals()
pcd.orient_normals_consistent_tangent_plane(7)

# mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=10, width=0, scale=50, linear_fit=False)[0]

mesh=o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd,7)

mesh.remove_degenerate_triangles()
mesh.remove_non_manifold_edges()

mesh.remove_duplicated_vertices()
mesh.remove_duplicated_triangles()

# mask=np.array([np.isin(Tri[:,0], label0).T, np.isin(Tri[:,1], label0).T,np.isin(Tri[:,2], label0).T]).T
# mask0=[]

# for i in range(0,len(mask)):
#     if mask[i].all() == True:
#         mask0.append([True])
#     elif mask[i].any() != True or mask[i].all() == False:
#         mask0.append([False])
# mask0=np.array(mask0)


triangle_clusters, cluster_n_triangles, cluster_area = (
            mesh.cluster_connected_triangles())
triangle_clusters = np.asarray(triangle_clusters)
cluster_n_triangles = np.asarray(cluster_n_triangles)
cluster_area = np.asarray(cluster_area)

mesh_0 = copy.deepcopy(mesh)
mesh_1= copy.deepcopy(mesh)
# triangles_to_remove = triangle_clusters>0
# mesh_0.remove_triangles_by_mask(mask0)
# mesh_0.remove_unreferenced_vertices()


o3d.visualization.draw_geometries([pcd,mesh_0], mesh_show_back_face=True)


# mesh = o3d.geometry.TriangleMesh.create_sphere(radius=0.61)
mesh = mesh.filter_smooth_taubin(number_of_iterations=5,lambda_filter=0.5)

mesh.paint_uniform_color([0.5, 0.5, 0.5])
#mesh.merge_close_vertices(1)
# mesh.remove_degenerate_triangles()
# mesh.remove_non_manifold_edges()

# mesh.remove_duplicated_vertices()
# mesh.remove_duplicated_triangles()
mesh.compute_vertex_normals()
mesh.compute_triangle_normals(normalized=True)
mesh.orient_triangles()


half_edge_mesh = o3d.geometry.HalfEdgeTriangleMesh.create_from_triangle_mesh(mesh)

o3d.visualization.draw_geometries([pcd,mesh], mesh_show_back_face=True)
if mesh.is_watertight()==True:
    print('Volume is :')
    print(mesh.get_volume())    
print('Surface Area is :')
print(mesh.get_surface_area())
mesh.is_watertight()

HalfEdge=np.asarray(half_edge_mesh.half_edges)
HalfEdge_normal=np.asarray(half_edge_mesh.triangle_normals)
HalfEdge_vertices=np.asarray(half_edge_mesh.vertices)
UsedEdges=[]
TurnAng=[]
MeanWidth=0
for i in range (0, len(HalfEdge)):
    
    if (i not in UsedEdges):
        CurrentEdge_Index=i
        TwinEdge_Index=HalfEdge[i].twin
        Current_Triangle_Index=HalfEdge[i].triangle_index
        Twin_Triangle_index= HalfEdge[TwinEdge_Index].triangle_index
        Current_Triangle_normal=HalfEdge_normal[Current_Triangle_Index]
        Twin_Triangle_normal=HalfEdge_normal[Twin_Triangle_index]
        DotProduct=np.dot(Current_Triangle_normal,Twin_Triangle_normal)
        
        if DotProduct > 1:
            DotProduct=1
        if DotProduct < -1:
            DotProduct=-1
            
        TurningAngle=np.arccos(DotProduct)
        Edge_vertex0=HalfEdge[i].vertex_indices[0]
        Edge_vertex1=HalfEdge[i].vertex_indices[1]
        Edge_Length= np.linalg.norm( HalfEdge_vertices[Edge_vertex0]-HalfEdge_vertices[Edge_vertex1])
        # print(Edge_Length)
        Edge_Length_Times_TurningAngle=Edge_Length*TurningAngle
        MeanWidth= MeanWidth + Edge_Length_Times_TurningAngle
        UsedEdges.append(TwinEdge_Index)
        TurnAng.append(TurningAngle)

MeanWidth= (1/(2*np.pi))*MeanWidth
print(MeanWidth)

tri_mesh = trimesh.Trimesh(np.asarray(mesh.vertices), np.asarray(mesh.triangles),
                          vertex_normals=np.asarray(mesh.vertex_normals))
tri_mesh.process
print(tri_mesh.integral_mean_curvature)
print(np.sqrt(4.0*np.pi*mesh.get_surface_area()))
# tri_mesh_0.show(viewer='gl')

# pcd_mesh = o3d.geometry.PointCloud()
# pcd_mesh.points = o3d.utility.Vector3dVector(np.asarray(mesh.vertices))
# pcd_mesh.estimate_normals()
# pcd_mesh.orient_normals_consistent_tangent_plane(10)
# o3d.visualization.draw_geometries([pcd_mesh])

'''

def xyz_to_pcd(xyz):
    pcd=o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(xyz)
    
    return(pcd)

def UnWrapCluster(pcd,EdgeLength):
    # NumberOfPoints=len(pcd.points)
    puc=copy.deepcopy(pcd)
    labels = np.array(puc.cluster_dbscan(eps=5, min_points=1))
    max_label = labels.max()
    id_big=np.where(labels == 0)[0]
    pcd_big = puc.select_by_index(id_big)
    d=[]
    for i in range(0,max_label+1):
        id_i=np.where(labels == i)[0]
        pcd_i = puc.select_by_index(id_i)
        diff_center=np.rint((pcd_big.get_center()-pcd_i.get_center())/EdgeLength)
        d.append(diff_center)
    
    pu_a=np.asarray(puc.points)
    pu_ac=copy.deepcopy(pu_a)
    for i in range(0,len(pu_ac)):
        
        l=labels[i]
        pu_ac[i]=pu_ac[i]+d[l]*EdgeLength
    return(xyz_to_pcd(pu_ac))
m=[]
for i in range(0,len(GrainBoundary)):
    
    if 2 in GrainBoundary[i][0] and len(GrainBoundary[i][1:])>50:
        k=i
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(GrainBoundary[k][1:])
        pcd=UnWrapCluster(pcd,EdgeLength)
        pcd.estimate_normals()
        pcd.orient_normals_consistent_tangent_plane(10)
        mesh, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
                pcd, depth=5)
        
        densities = np.asarray(densities)
        max_d=np.max(densities)
        density_colors = plt.get_cmap('plasma')(
            (densities - densities.min()) / (densities.max() - densities.min()))
        density_colors = density_colors[:, :3]
        density_mesh = o3d.geometry.TriangleMesh()
        density_mesh.vertices = mesh.vertices
        density_mesh.triangles = mesh.triangles
        density_mesh.triangle_normals = mesh.triangle_normals
        density_mesh.vertex_colors = o3d.utility.Vector3dVector(density_colors)
        # o3d.visualization.draw_geometries([density_mesh,pcd])
        
        # vertices_to_remove = densities < np.quantile(densities, 0.1)
        vertices_to_remove = densities < np.max(densities)-0.9
        density_mesh.remove_vertices_by_mask(vertices_to_remove)
        density_mesh.remove_degenerate_triangles()
        density_mesh.remove_duplicated_triangles()
        density_mesh.remove_duplicated_vertices()
        density_mesh.remove_non_manifold_edges()
        # o3d.visualization.draw_geometries([density_mesh,pcd],mesh_show_back_face=True)
        
        mesh=copy.deepcopy(density_mesh)
        
        triangle_clusters, cluster_n_triangles, cluster_area = (
               mesh.cluster_connected_triangles())
        
        triangle_clusters = np.asarray(triangle_clusters)
        cluster_n_triangles = np.asarray(cluster_n_triangles)
        cluster_area = np.asarray(cluster_area)
        
        mesh_0 = copy.deepcopy(mesh)
        max_n_triangle=np.max(cluster_n_triangles)
        triangles_to_remove = cluster_n_triangles[triangle_clusters] < max_n_triangle
        mesh_0.remove_triangles_by_mask(triangles_to_remove)
        mesh_0 = mesh_0.filter_smooth_taubin(number_of_iterations=15,lambda_filter=0.5)
        # o3d.visualization.draw_geometries([mesh_0])
        m.append(mesh_0)
o3d.visualization.draw_geometries(m,mesh_show_back_face=True)


k=27
pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(GrainBoundary[k][1:])
pcd=UnWrapCluster(pcd,EdgeLength)
pcd.estimate_normals()C
pcd.orient_normals_consistent_tangent_plane(10)

distances = pcd.compute_nearest_neighbor_distance()
avg_dist = np.mean(distances)
radius = 3 * avg_dist
bpa_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(pcd,o3d.utility.DoubleVector([radius, radius * 2]))
bpa_mesh.remove_degenerate_triangles()
bpa_mesh.remove_duplicated_triangles()
bpa_mesh.remove_duplicated_vertices()
bpa_mesh.remove_non_manifold_edges()
o3d.visualization.draw_geometries([bpa_mesh,pcd],mesh_show_back_face=True)

m1=[]
for i in range(0,len(GrainBoundary)):
    
    if 2 in GrainBoundary[i][0] and len(GrainBoundary[i][1:])>10:
        k=i
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(GrainBoundary[k][1:])
        pcd=UnWrapCluster(pcd,EdgeLength)
        mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, 4)
        mesh.paint_uniform_color([random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)])
        m1.append(mesh)
        
o3d.visualization.draw_geometries(m1, mesh_show_back_face=True)

gb0=copy.deepcopy(m[5])

tri_mesh = trimesh.Trimesh(np.asarray(gb0.vertices), np.asarray(gb0.triangles),
                          vertex_normals=np.asarray(gb0.vertex_normals))

tri_mesh.fix_normals()
thick=trimesh.proximity.thickness(tri_mesh, tri_mesh.vertices, exterior=False, normals=None, method='ray')
vn=tri_mesh.vertex_normals
t=np.expand_dims(thick, axis=-1)
tri_mesh.vertices=tri_mesh.vertices-1*t*vn
tri_mesh=tri_mesh.process()
tri_mesh.show(viewer='gl')