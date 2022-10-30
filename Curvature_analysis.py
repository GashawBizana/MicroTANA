import numpy as np
import open3d as o3d
import trimesh
import copy
import time
import pyvista as pv

import time
# def GrainAtoms():
#     GrainAtom=[[]for i in max(GrainId)]
#     for i in range(0,len(x)):
        
#         GrainAtom.appendGrainId[i]==i
  
        

def ConstructPeriodicLocation(P, EdgeLength):
    
    p0=P[0]
    p1=P[1]
    p2=P[2]
    p0m=P[0]-EdgeLength[0]
    p1m=P[1]-EdgeLength[1]
    p2m=P[2]-EdgeLength[2]
    p0m=P[0]-EdgeLength[0]
    p1m=P[1]-EdgeLength[1]
    p2m=P[2]-EdgeLength[2]
    p0p=P[0]+EdgeLength[0]
    p1p=P[1]+EdgeLength[1]
    p2p=P[2]+EdgeLength[2]
    p0_l=[p0,p0m,p0p]
    p1_l=[p1,p1m,p1p]
    p2_l=[p2,p2m,p2p]
    
    P_periodic=[]
    
    for i in range(0,3):
        
        for j in range(0,3):
            
            for k in range(0,3):
                
                P_periodic.append([p0_l[i],p1_l[j],p2_l[k]])
    return(P_periodic)

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
#     return(pcd_all)

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
    pcdout=xyz_to_pcd(pu_ac)
    return(pcdout)

def xyz_to_pcd(xyz):
    pcd=o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(xyz)
    
    return(pcd)



def GrainMesh (xyz,EdgeLength):
    pcdWrapped = xyz_to_pcd(xyz)
    pcd=UnWrapCluster(pcdWrapped,EdgeLength)
    mesh=o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd,4)
    mesh.remove_degenerate_triangles()
    mesh.remove_non_manifold_edges()
    mesh.remove_duplicated_vertices()
    mesh.remove_duplicated_triangles()
    mesh.orient_triangles()
    mesh_smoothed=copy.deepcopy(mesh)#
    # mesh_smoothed=o3d.geometry.TriangleMesh.filter_smooth_taubin(mesh_smoothed, number_of_iterations=1,lambda_filter=0.5,mu=-0.526315789)
    
    return(mesh,mesh_smoothed)

def GrainBondaryMesh1(mesh_grain,mesh_smoothed_i,X_GB,EdgeLength):
    mesh_grain_temp=copy.deepcopy(mesh_grain)
    mesh_grain_smoothed_temp=copy.deepcopy(mesh_smoothed_i)
    triangle_vertices=np.asarray(mesh_grain_temp.vertices)
    verteices_to_remove=[]
    X_GB=np.array(X_GB)
    for i in range(0,len(triangle_vertices)):
        # triangel_vertices_periodic_location=ConstructPeriodicLocation(triangle_vertices[i], EdgeLength)
        # # triangel_vertices_periodic_location=triangel_vertices_periodic_location[0]
        
        # remove_condition=[]
        
        # for index in range (0,len(triangel_vertices_periodic_location)):
            
        #     con= triangel_vertices_periodic_location[index] in X_GB#((X_GB == triangel_vertices_periodic_location[index]).all(1).any())
        #     #con=np.isin(triangel_vertices_periodic_location[index],X_GB)
        #     remove_condition.append(con)
        
        
        # if any(remove_condition)==False:
        if ((np.abs(triangle_vertices[i]-X_GB))%EdgeLength).all() > (20/min(EdgeLength)):
            verteices_to_remove.append(i)
            
    # mesh_grain_temp.remove_vertices_by_index(verteices_to_remove)
    mesh_grain_smoothed_temp.remove_vertices_by_index(verteices_to_remove)
    mesh_grain_smoothed_temp.remove_unreferenced_vertices()
    
    mesh_GB_smoothed=mesh_grain_smoothed_temp
    

    mesh_GB_smoothed.compute_vertex_normals()
    
    return(mesh_GB_smoothed)

def GrainBondaryMesh2(mesh_grain,mesh_smoothed_i,X_GB,EdgeLength):
    mesh_grain_temp=copy.deepcopy(mesh_grain)
    mesh_grain_smoothed_temp=copy.deepcopy(mesh_smoothed_i)
    triangle_vertices=np.asarray(mesh_grain_temp.vertices)
    triangles=np.asarray(mesh_grain_temp.triangles)
    triangles_to_remove=[]
    X_GB=np.array(X_GB)
    for i in range(0,len(triangles)):
        p0,p1,p2=triangle_vertices[triangles[i][0]],triangle_vertices[triangles[i][1]],triangle_vertices[triangles[i][2]]
        Test=[((np.abs(p0-X_GB))%EdgeLength).all() > 4/min(EdgeLength), ((np.abs(p1-X_GB))%EdgeLength).all() > 4/min(EdgeLength),((np.abs(p2-X_GB))%EdgeLength).all() > 4/min(EdgeLength)]
        if Test.count(False)<1:
        # triangel_vertices_periodic_location=ConstructPeriodicLocation(triangle_vertices[i], EdgeLength)
        # # triangel_vertices_periodic_location=triangel_vertices_periodic_location[0]
        
        # remove_condition=[]
        
        # for index in range (0,len(triangel_vertices_periodic_location)):
            
        #     con= triangel_vertices_periodic_location[index] in X_GB#((X_GB == triangel_vertices_periodic_location[index]).all(1).any())
        #     #con=np.isin(triangel_vertices_periodic_location[index],X_GB)
        #     remove_condition.append(con)
        
        
        # if any(remove_condition)==False:
        # if ((np.abs(triangle_vertices[i]-X_GB))%EdgeLength).all() > (20/min(EdgeLength)):
            triangles_to_remove.append(i)
            
    # mesh_grain_temp.remove_vertices_by_index(verteices_to_remove)
    mesh_grain_smoothed_temp.remove_triangles_by_index(triangles_to_remove)
    mesh_grain_smoothed_temp.remove_unreferenced_vertices()
    
    mesh_GB_smoothed=mesh_grain_smoothed_temp
    

    mesh_GB_smoothed.compute_vertex_normals()
    
    return(mesh_GB_smoothed)


def GrainBondaryMesh(mesh_grain,mesh_smoothed_i,X_GB,EdgeLength,InnerAtom_not_gb):
    
    # 
    def Test_vertice(diff,EdgeLength):
        is_vertice_part_of_mesh=np.any(np.all(np.logical_or( np.logical_and(np.less_equal(0,diff), np.less_equal(diff,0)) , np.logical_and(np.less_equal(EdgeLength-0,diff), np.less_equal(diff,EdgeLength+0)) ), axis=1))
        return(is_vertice_part_of_mesh)
    
    # def test(dx,EdgeLength):
    #     test=[ ((0 <= i <=4) or (EdgeLength-4 <= i <= EdgeLength + 4)) for i in dx ]
        
    #     return(test)
        
    
    mesh_grain_temp=copy.deepcopy(mesh_grain)
    mesh_grain_smoothed_temp=copy.deepcopy(mesh_smoothed_i)
    triangle_vertices=np.asarray(mesh_grain_temp.vertices)
    triangles=np.asarray(mesh_grain_temp.triangles)
    triangles_to_remove=[]
    X_GB=np.array(X_GB)
    InnerAtom_not_gb=np.array(InnerAtom_not_gb)
    for i in range(0,len(triangles)):
        p0,p1,p2=triangle_vertices[triangles[i][0]],triangle_vertices[triangles[i][1]],triangle_vertices[triangles[i][2]]
        # p00,p01,p02=p0[:,0],p0[:,1],p0[:,2]
        # p10,p11,p12=p1[:,0],p1[:,1],p1[:,2]
        # p20,p21,p22=p2[:,0],p2[:,1],p2[:,2]
        
        # dx0,dx1,dx2=np.abs(p0-X_GB),np.abs(p1-X_GB),np.abs(p2-X_GB)
        
        # dx00,dx01,dx02=dx0[:,0],dx0[:,1],dx0[:,2]
        # dx10,dx11,dx12=dx1[:,0],dx1[:,1],dx1[:,2]
        # dx20,dx21,dx22=dx2[:,0],dx2[:,1],dx2[:,2]
        
        # dx00_T,dx01_T,dx02_T=test(dx00,EdgeLength[0]),test(dx01,EdgeLength[1]),test(dx02,EdgeLength[2])
        
        # dx0_T=any([all([dx00_T[idx],dx01_T[idx],dx02_T[idx]]) for idx,j in enumerate(dx00_T)])
        
        # dx10_T,dx11_T,dx12_T=test(dx10,EdgeLength[0]),test(dx11,EdgeLength[1]),test(dx12,EdgeLength[2])
        
        # dx1_T=any([all([dx10_T[idx],dx11_T[idx],dx12_T[idx]]) for idx,j in enumerate(dx10_T)])
        
        # dx20_T,dx21_T,dx22_T=test(dx20,EdgeLength[0]),test(dx21,EdgeLength[1]),test(dx22,EdgeLength[2])
        
        # dx2_T=any([all([dx20_T[idx],dx21_T[idx],dx22_T[idx]]) for idx,j in enumerate(dx20_T)])
        
        
        Test=[Test_vertice(np.abs(p0-X_GB),EdgeLength),Test_vertice(np.abs(p1-X_GB),EdgeLength),Test_vertice(np.abs(p2-X_GB),EdgeLength)]
        Test2=[Test_vertice(np.abs(p0-InnerAtom_not_gb),EdgeLength),Test_vertice(np.abs(p1-InnerAtom_not_gb),EdgeLength),Test_vertice(np.abs(p2-InnerAtom_not_gb),EdgeLength)]
        if Test.count(True) <=1: 
            if Test2.count(True) <2 :
        # triangel_vertices_periodic_location=ConstructPeriodicLocation(triangle_vertices[i], EdgeLength)
        # # triangel_vertices_periodic_location=triangel_vertices_periodic_location[0]
        
        # remove_condition=[]
        
        # for index in range (0,len(triangel_vertices_periodic_location)):
            
        #     con= triangel_vertices_periodic_location[index] in X_GB#((X_GB == triangel_vertices_periodic_location[index]).all(1).any())
        #     #con=np.isin(triangel_vertices_periodic_location[index],X_GB)
        #     remove_condition.append(con)
        
        
        # if any(remove_condition)==False:
        # if ((np.abs(triangle_vertices[i]-X_GB))%EdgeLength).all() > (20/min(EdgeLength)):
                triangles_to_remove.append(i)
            
    # mesh_grain_temp.remove_vertices_by_index(verteices_to_remove)
    mesh_grain_smoothed_temp.remove_triangles_by_index(triangles_to_remove)
    mesh_grain_smoothed_temp.remove_unreferenced_vertices()
    
    mesh_GB_smoothed=mesh_grain_smoothed_temp
    

    mesh_GB_smoothed.compute_vertex_normals()
    
    return(mesh_GB_smoothed)

def GrainBoundaryCurvatureAndSurfaceArea(mesh_grain_boundary):
    tri_mesh_GB = trimesh.Trimesh(np.asarray(mesh_grain_boundary.vertices), np.asarray(mesh_grain_boundary.triangles),
                              vertex_normals=np.asarray(mesh_grain_boundary.vertex_normals), triangle_normals=np.asarray(mesh_grain_boundary.vertex_normals))
    tri_mesh_GB.process(validate=True, merge_tex=None, merge_norm=None)
    
    trimesh.repair.fill_holes(tri_mesh_GB)
    trimesh.repair.fix_normals(tri_mesh_GB, multibody=False)
    trimesh.repair.fix_inversion(tri_mesh_GB, multibody=False)
    
    total_mean_curvature=tri_mesh_GB.integral_mean_curvature
    surface_area=tri_mesh_GB.area
    average_gb_mean_curvature=total_mean_curvature/surface_area
    return(total_mean_curvature,surface_area,average_gb_mean_curvature,tri_mesh_GB)



def GrainBoundaryProperties(grain_boundary_list,grain_inner_atom_list,EdgeLength,X,GrainId,ModifiedGrainIdListSorted):
    number_of_grain_boundaries=len(grain_boundary_list)
    number_of_grains=len(grain_inner_atom_list)
    grain_boundary_property=[]
    grain_property=[]
    
    for i in range(0,number_of_grains):
        InnerAtom_not_gb=grain_inner_atom_list[i][1:]
        mid=[(ModifiedGrainIdListSorted[idx][2]).tolist() for idx,j in enumerate(X) if GrainId[idx] ==i]
        gid=[GrainId[idx] for idx,j in enumerate(X) if GrainId[idx] ==i]
        InnerAtom_i=[X[Idx].tolist() for Idx in range(len(X)) if GrainId[Idx]==i]
        mesh_i,mesh_smoothed_i=GrainMesh(InnerAtom_i,EdgeLength)
        mesh_i_Trimesh=trimesh.Trimesh(vertices=np.asarray(mesh_i.vertices), faces=np.asarray(mesh_i.triangles))
        
        grain_property.append([i,mesh_i_Trimesh])
        for j in range (0, number_of_grain_boundaries):
            if i in grain_boundary_list[j][0]: #(i== grain_boundary_list[j][0][0] or i==grain_boundary_list[j][0][1]):
                grain_boundary_ij= grain_boundary_list[j][1:]
                grain_boundary_mesh_ij=GrainBondaryMesh(mesh_i,mesh_smoothed_i,grain_boundary_ij,EdgeLength,InnerAtom_not_gb)
               
                
                grain_boundary_mesh_Trimesh=trimesh.Trimesh(vertices=np.asarray(grain_boundary_mesh_ij.vertices), faces=np.asarray(grain_boundary_mesh_ij.triangles))
                

                grain_boundary_property.append([[i,[k for k in grain_boundary_list[j][0] if k!=i][0]],grain_boundary_mesh_Trimesh])
                
    return(grain_boundary_property,grain_property)



def GrainBoundaryPropertyExtended(GrainBoundary,grain_boundary_property_extended):
    
    grain_boundary_property=[]
    
    for i in range (0,len(GrainBoundary)):
        gb_Id=GrainBoundary[i][0]
        
        gb_meshes_extended=[j[1] for j in grain_boundary_property_extended if ((gb_Id[0] in j[0]) and (gb_Id[1] in j[0]))]
        
        grain_boundary_property.append([gb_Id,gb_meshes_extended])
        
    return(grain_boundary_property)
        
        



    
# Test=[np.all(np.greater_equal(np.abs(m.triangles[i][0]-gb0)%EdgeLength,4/min(EdgeLength))),np.all(np.greater_equal(np.abs(m.triangles[i][1]-gb0)%EdgeLength,4/min(EdgeLength))),np.all(np.greater_equal(np.abs(m.triangles[i][2]-gb0)%EdgeLength,4/min(EdgeLength)))]


def test_for_clossnes(diff,tol,EdgeLength):
    Test0=[np.logical_or(np.logical_and(np.greater_equal(diff[:,0],-tol), np.less_equal(diff[:,0],tol)),np.logical_and(np.greater_equal(diff[:,0],EdgeLength[0]-tol),np.less_equal(diff[:,0],EdgeLength[0]+tol))),
           np.logical_or(np.logical_and(np.greater_equal(diff[:,1],-tol), np.less_equal(diff[:,1],tol)),np.logical_and(np.greater_equal(diff[:,1],EdgeLength[1]-tol),np.less_equal(diff[:,1],EdgeLength[1]+tol))),
           np.logical_or(np.logical_and(np.greater_equal(diff[:,2],-tol), np.less_equal(diff[:,2],tol)),np.logical_and(np.greater_equal(diff[:,2],EdgeLength[2]-tol),np.less_equal(diff[:,2],EdgeLength[2]+tol)))]
    

    return(np.any(np.all(Test0,axis=0)))
# Xgb=np.array(grain_boundary_list[j][1:])
# m=copy.deepcopy(m2)
# tol=1.5*1.43
# T=m.triangles
# mask=[]
# mask=[]
# triangle=[]
# for i in range(0,len(m.triangles)):
#     v0,v1,v2=T[i][0],T[i][1],T[i][2]
#     diff0,diff1,diff2=np.abs(v0-Xgb)%EdgeLength,np.abs(v1-Xgb)%EdgeLength,np.abs(v2-Xgb)%EdgeLength
#     Test=[test_for_clossnes(diff0,tol,EdgeLength),test_for_clossnes(diff1,tol,EdgeLength),test_for_clossnes(diff2,tol,EdgeLength)]
#     if Test.count(True) >=3:
#         mask.append(True)
#         triangle.append(T[i])
        
#     else :
#         mask.append(False)
# m.update_faces(mask)
# m.show(viewer='gl')


# trimesh.proximity.closest_point_naive(m2, xgb)
