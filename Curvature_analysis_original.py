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
        print(i)
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
    mesh=o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd,3.5)
    mesh.remove_degenerate_triangles()
    mesh.remove_non_manifold_edges()
    mesh.remove_duplicated_vertices()
    mesh.remove_duplicated_triangles()
    mesh.orient_triangles()
    mesh_smoothed=copy.deepcopy(mesh)#o3d.geometry.TriangleMesh.filter_smooth_taubin(mesh, number_of_iterations=0,lambda_filter=0.5,mu=-0.526315789)
    
    return(mesh,mesh_smoothed)

def GrainBondaryMesh(mesh_grain,mesh_smoothed_i,X_GB,EdgeLength):
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
        if ((np.abs(triangle_vertices[i]-X_GB))%EdgeLength).all() > 0.00001:
            verteices_to_remove.append(i)
            
    # mesh_grain_temp.remove_vertices_by_index(verteices_to_remove)
    mesh_grain_smoothed_temp.remove_vertices_by_index(verteices_to_remove)
    mesh_grain_smoothed_temp.remove_non_manifold_edges()
    mesh_grain_smoothed_temp.remove_unreferenced_vertices()
    #mesh_GB=o3d.geometry.TriangleMesh.filter_smooth_taubin(mesh_grain_temp, number_of_iterations=0,lambda_filter=0.5, mu=-0.526315789)
    # mesh_GB=mesh_grain_temp
    mesh_GB_smoothed=mesh_grain_smoothed_temp
    
    # mesh_GB.compute_vertex_normals()
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



def GrainBoundaryProperties(grain_boundary_list,grain_inner_atom_list,EdgeLength,X,GrainId):
    number_of_grain_boundaries=len(grain_boundary_list)
    number_of_grains=len(grain_inner_atom_list)
    grain_boundary_property=[]
    grain_property=[]
    
    for i in range(0,number_of_grains):
        InnerAtom_i=np.array(grain_inner_atom_list[i][1:])
        InnerAtom_i=[X[Idx].tolist() for Idx in range(len(X)) if GrainId[Idx]==i]
        mesh_i,mesh_smoothed_i=GrainMesh(InnerAtom_i,EdgeLength)
        mesh_smoothed_i_Trimesh=trimesh.Trimesh(np.asarray(mesh_smoothed_i.vertices), np.asarray(mesh_smoothed_i.triangles),
                                  vertex_normals=np.asarray(mesh_smoothed_i.vertex_normals))
        grain_property.append([i,mesh_smoothed_i_Trimesh])
        for j in range (0, number_of_grain_boundaries):
            if i in grain_boundary_list[j][0]: #(i== grain_boundary_list[j][0][0] or i==grain_boundary_list[j][0][1]):
                grain_boundary_ij= grain_boundary_list[j][1:]
                grain_boundary_mesh_ij=GrainBondaryMesh(mesh_i,mesh_smoothed_i,grain_boundary_ij,EdgeLength)
               
                
                if grain_boundary_mesh_ij.has_triangles():
                    
                    grain_boundary_mesh_ij_vol,grain_boundary_mesh_smoothed_ij_vol=GrainMesh (grain_boundary_ij,EdgeLength)
                    grain_boundary_mesh_smoothed_ij_vol_Trimesh=trimesh.Trimesh(np.asarray(grain_boundary_mesh_smoothed_ij_vol.vertices), np.asarray(grain_boundary_mesh_smoothed_ij_vol.triangles),
                                              vertex_normals=np.asarray(grain_boundary_mesh_smoothed_ij_vol.vertex_normals))
                    
                    total_mean_curvature_ij,surface_area_ij,average_gb_mean_curvature_ij,tri_mesh_GB_ij=GrainBoundaryCurvatureAndSurfaceArea(grain_boundary_mesh_ij)
                    # grain_boundary_property.append([grain_boundary_list[j][0],total_mean_curvature_ij,surface_area_ij,grain_boundary_mesh_ij])
                    
                    
                    # grain_boundary_property.append([[i,[k for k in grain_boundary_list[j][0] if k!=i][0]],total_mean_curvature_ij,surface_area_ij,np.array(grain_boundary_mesh_ij.vertices)])
                    grain_boundary_property.append([[i,[k for k in grain_boundary_list[j][0] if k!=i][0]],total_mean_curvature_ij,surface_area_ij,average_gb_mean_curvature_ij,tri_mesh_GB_ij,grain_boundary_mesh_smoothed_ij_vol_Trimesh])
                    
    return(grain_boundary_property,grain_property)



def GrainBoundaryPropertyExtended(GrainBoundary,grain_boundary_property_extended):
    
    grain_boundary_property=[]
    
    for i in range (0,len(GrainBoundary)):
        gb_Id=GrainBoundary[i][0]
        
        gb_meshes_extended=[j[1] for j in grain_boundary_property_extended if ((gb_Id[0] in j[0]) and (gb_Id[1] in j[0]))]
        
    grain_boundary_property.append()
        
        

# def GrainProperties(grain_inner_atom_list):
    
    

# def QuantityCalculation(mesh_grain, grain_boundary_list,grain_id, EdgeLength):
#     grain_boundaries_i= [grain_boundary[1:] for grain_boundary in grain_boundary_list if grain_id in grain_boundary[0]]
#     grains_forming_grain_boundary=[grain_boundary[0] for grain_boundary in grain_boundary_list if grain_id in grain_boundary[0]]
#     grain_boundaries_mesh_i= [GrainBondaryMesh(mesh_grain,grain_boundary_ij,EdgeLength) for grain_boundary_ij in grain_boundaries_i ]
#     total_mean_curvature_and_surface_area_i= [GrainBoundaryCurvatureAndSurfaceArea(grain_boundary_mesh_ij) for grain_boundary_mesh_ij in grain_boundaries_mesh_i ]
    
#     return(grains_forming_grain_boundary,total_mean_curvature_and_surface_area_i,grain_boundaries_mesh_i)
# def GrainBoundaryProperties(grain_boundary_list,grain_inner_atom_list,EdgeLength):
#     number_of_grain_boundaries=len(grain_boundary_list)
#     number_of_grains=len(grain_inner_atom_list)
#     grain_boundary_property=[]
#     for i in range(0,1):
#         InnerAtom_i=np.array(grain_inner_atom_list[i][1:])
        
#         mesh_i,mesh_smoothed_i=GrainMesh(InnerAtom_i,EdgeLength)
        
#         grain_boundaries_i= [grain_boundary[1:] for grain_boundary in grain_boundary_list if i in grain_boundary[0]]
#         grain_boundaries_mesh_i= [GrainBondaryMesh(mesh_i,grain_boundary_ij,EdgeLength) for grain_boundary_ij in grain_boundaries_i ]
#         total_mean_curvature_and_surface_area_i= [GrainBoundaryCurvatureAndSurfaceArea(grain_boundary_mesh_ij) for grain_boundary_mesh_ij in grain_boundaries_mesh_i ]
#         for j in range (0, number_of_grain_boundaries):
#             if (i== grain_boundary_list[j][0][0] or i==grain_boundary_list[j][0][1]):
#                 grain_boundary_ij= grain_boundary_list[j][1:]
#                 grain_boundary_mesh_ij=GrainBondaryMesh(mesh_i,grain_boundary_ij,EdgeLength)
#                 total_mean_curvature_ij,surface_area_ij=GrainBoundaryCurvatureAndSurfaceArea(grain_boundary_mesh_ij)
        # grain_boundary_property.append([[i,j],[total_mean_curvature_ij],[surface_area_ij],[grain_boundary_mesh_ij]])
#     return(grain_boundary_property)

def GrainBoundaryDisplacement_using_mesh(grain_boundary_n_mesh, grain_boundary_m_mesh,EdgeLength):
    
    center_of_mass_n=grain_boundary_n_mesh.center_mass
    center_of_mass_m=grain_boundary_m_mesh.center_mass
    
    diff_center=np.rint((center_of_mass_n-center_of_mass_m)/EdgeLength)    
    periodic_translation_factor=diff_center*EdgeLength
    
    transformation_matrix,cost_rmsd=trimesh.registration.mesh_other(grain_boundary_n_mesh,grain_boundary_m_mesh, samples=2000,scale=False,icp_first=50,icp_final=250)
    
    translation_based_on_cm=center_of_mass_n-center_of_mass_m+ periodic_translation_factor
                
    return(transformation_matrix,cost_rmsd,translation_based_on_cm)

def DistanceBetweenSurfaces(surf1,surf2):
    p = pv.Plotter()
    p.add_mesh(surf1, smooth_shading=True)
    p.add_mesh(surf2, smooth_shading=True)
    p.show_grid()
    p.show()


def GrainIdMapper(grain_mapping_data, grain_boundary_forming_grains_n,time_step_n):

    n0=grain_boundary_forming_grains_n[0]
    n1=grain_boundary_forming_grains_n[1]
    grain_mapping_data_time_step_n=grain_mapping_data[time_step_n]
    index_n0=np.where(grain_mapping_data_time_step_n==n0)
    index_n1=np.where(grain_mapping_data_time_step_n==n1)
    grain_mapping_data_time_step_m=grain_mapping_data[time_step_n-1]
    m0=grain_mapping_data_time_step_m[index_n0]
    m1=grain_mapping_data_time_step_m[index_n1]

    return([int(m0[0]),int(m1[0])])
 
def GrainBoundaryVelocity(grain_boundary_property_n, grain_boundary_property_m,EdgeLength):
    
    grain_boundary_velocity_curvature_area=[]
    for i in range (0,len(grain_boundary_property_n)):
        
        grain_boundary_forming_grains_n=grain_boundary_property_n[i][0]
        grain_boundary_forming_grains_m=GrainIdMapper(grain_boundary_forming_grains_n)
        
        for j in range(0,len(grain_boundary_property_m)):
            
            if grain_boundary_property_m[j][0]== grain_boundary_forming_grains_m:
                
                transformation_matrix,cost_rmsd,translation_based_on_cm=GrainBoundaryDisplacement_using_mesh(grain_boundary_property_n[i][3], grain_boundary_property_m[j][3],EdgeLength)
                
                grain_boundary_velocity_curvature_area.append([[grain_boundary_forming_grains_n],[transformation_matrix],[translation_based_on_cm]])
    return(grain_boundary_velocity_curvature_area)





# t = time.time()
# grain_boundary_property_current=GrainBoundaryProperties(GrainBoundary,InnerAtom,EdgeLength)
# print ("Curvature analysis for one time step took took: ", (time.time()-t))
# grain_boundary_property_next=GrainBoundaryProperties(GrainBoundary_1,InnerAtom_1,EdgeLength_1)

    
    