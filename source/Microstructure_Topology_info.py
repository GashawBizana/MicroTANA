import numpy as np
import trimesh
import pyvista as pv
import pymeshfix as mf
from pymeshfix._meshfix import PyTMesh

from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
from ovito.vis import *
from ovito.qt_compat import QtCore

from sklearn.cluster import DBSCAN
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

def Grain_mesh_ovito(GrainId,X,EdgeLength,grainId):
    
    data = DataCollection()
    cell = SimulationCell(pbc = (False, False, False))
    cell[:,0] = (EdgeLength[0],0,0)
    cell[:,1] = (0,EdgeLength[1],0)
    cell[:,2] = (0,0,EdgeLength[2])
    data.objects.append(cell)
    pos_u=[X[idx] for idx,i in enumerate(X) if GrainId[idx] ==grainId]
    pos=UnWrapCluster(np.array(pos_u),EdgeLength)

    gid=[GrainId[idx] for idx,i in enumerate(X) if GrainId[idx] ==grainId]
    particles = Particles()
    particles.create_property('Position', data=pos)
    particles.create_property('grainId', data=gid)
    data.objects.append(particles)

    # Create a new Pipeline with a StaticSource as data source:
    pipeline = Pipeline(source = StaticSource(data = data))

    modifier_select=ExpressionSelectionModifier(expression = f"grainId=={grainId}")

    pipeline.modifiers.append(modifier_select)
    modifier_surface=ConstructSurfaceModifier(
        method = ConstructSurfaceModifier.Method.GaussianDensity, grid_resolution=100,
        radius_scaling = 1,
        isolevel = 0.9,only_selected=True)
    pipeline.modifiers.append(modifier_surface)
    # pipeline.modifiers.append(ConstructSurfaceModifier(radius = 7,only_selected=True))
    data1=pipeline.compute()
    mesh = data1.surfaces['surface']
    tri_mesh=trimesh.Trimesh(vertices=mesh.get_vertices(),faces=mesh.get_faces())
    return(tri_mesh)


def misorientation_quaternions(timestep,grain_boundary_forming_step_i, q0,q1,q2,q3,grain_mapping_data):
    
    n0_before_mapping=grain_boundary_forming_step_i[0] # before mapping
    n1_before_mapping=grain_boundary_forming_step_i[1] # before mapping
    
    grain_mapping_data_time_step_n=grain_mapping_data[timestep]
    n0=np.where(grain_mapping_data_time_step_n==n0_before_mapping)
    n1=np.where(grain_mapping_data_time_step_n==n1_before_mapping)
    
    q0_data_time_step_n=q0[timestep]
    q1_data_time_step_n=q1[timestep]
    q2_data_time_step_n=q2[timestep]
    q3_data_time_step_n=q3[timestep]
    
    q0_n0=q0_data_time_step_n[n0]
    q1_n0=q1_data_time_step_n[n0]
    q2_n0=q2_data_time_step_n[n0]
    q3_n0=q3_data_time_step_n[n0]
    
    q0_n1=q0_data_time_step_n[n1]
    q1_n1=q1_data_time_step_n[n1]
    q2_n1=q2_data_time_step_n[n1]
    q3_n1=q3_data_time_step_n[n1]
    
    q1=np.array([q0_n0,q1_n0,q2_n0,q3_n0])
    q2=np.array([q0_n1,q1_n1,q2_n1,q3_n1])
    
    #misorientation_grain_boundary_forming_step_i= np.degrees( np.abs(    2* (np.arccos (q0_n0*q0_n1 + q1_n0*q1_n1 + q2_n0*q2_n1 + q3_n0*q3_n1 ))   )  )
    
    
    
    return(q1,q2)
    
def cubically_equivalent_quaternions (q1,q2):
    
    def q_mult(q1, q2):
        """
        
        
        Parameters
        ----------
        q1, q1: numpy ndarray(4)
            Quaternions. First item is the scalar part
        
        Returns
        -------
        q : numpy ndarray(4)
            Quaternion product
        """
        q = np.zeros(4)
        q[0] = q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2] + q1[3]*q2[3]
        q[1] = -q1[0]*q2[1] + q1[1]*q2[0] - q1[2]*q2[3] + q1[3]*q2[2]
        q[2] = -q1[0]*q2[2] + q1[1]*q2[3] + q1[2]*q2[0] - q1[3]*q2[1]
        q[3] = -q1[0]*q2[3] -q1[1]*q2[2] + q1[2]*q2[1] + q1[3]*q2[0]
        return (q)

    
    
    q=q_mult(q1, q2)
    
    a0,a1,a2,a3=q[0],q[1],q[2],q[3]
    
    cubic_sym=np.array([
    np.array([a0, a1, a2, a3]),
    
    np.array([-a1, a0, a3,-a2]),
    
    np.array([-a2,- a3, a0, a1]),
    
    np.array([-a3, a2, - a1, a0]),
    
    0.5*np.array([a0 - a1 - a2 - a3, a0 + a1 + a2 - a3, a0 - a1 + a2 + a3, a0 + a1 - a2 + a3]) ,
    
    0.5*np.array([a0 + a1 + a2 + a3,-a0 + a1 - a2 + a3, - a0 + a1 + a2 - a3, - a0 - a1 + a2 + a3]),
    
    0.5*np.array([a0 - a1 + a2 - a3, a0 + a1 + a2 + a3, -a0 - a1 + a2 + a3, a0 - a1 - a2 + a3]) ,
    
    0.5*np.array([a0 + a1 - a2 + a3,-a0 + a1 - a2 - a3, a0 + a1 + a2 - a3, - a0 + a1 + a2 + a3]) ,
    
    0.5*np.array([a0 + a1 - a2 - a3,-a0 + a1 + a2 - a3, a0 - a1 + a2 - a3, a0 + a1 + a2 + a3] ),
    
    0.5*np.array([a0 - a1 + a2 + a3, a0 + a1 - a2 + a3,-a0 + a1 + a2 + a3, - a0 - a1 - a2 + a3]),
    
    0.5*np.array([a0 + a1 + a2 - a3,-a0 + a1 + a2 + a3,-a0 - a1 + a2 - a3, a0 - a1 + a2 + a3]),
    
    0.5*np.array([a0 - a1 - a2 + a3, a0 + a1 - a2 - a3, a0 + a1 + a2 + a3, - a0 + a1 - a2 + a3]),
    
    (1/np.sqrt(2))*np.array([a0 - a1, a0 + a1, a2 + a3,-a2 + a3]) ,
    (1/np.sqrt(2))*np.array([a0 - a2, a1 - a3, a0 + a2, a1 + a3]) ,
    (1/np.sqrt(2))*np.array([a0 - a3, a1 + a2, -a1 + a2, a0 + a3]) ,
    (1/np.sqrt(2))*np.array([-a1 - a2, a0 - a3, a0 + a3, a1 - a2]),
    (1/np.sqrt(2))*np.array([-a2 - a3, a2 - a3, a0 - a1, a0 + a1]) ,
    (1/np.sqrt(2))*np.array([-a1 - a3, a0 + a2, -a1 + a3, a0 - a2]),
    (1/np.sqrt(2))*np.array([a0 + a1, -a0 + a1, a2 - a3, a2 + a3]),
    (1/np.sqrt(2))*np.array([a0 + a2, a1 + a3, -a0 + a2,- a1 + a3]),
    (1/np.sqrt(2))*np.array([a0 + a3, a1 - a2, a1 + a2, - a0 + a3]) ,
    (1/np.sqrt(2))*np.array([a1 - a2, - a0 - a3, a0 - a3, a1 + a2] ),
    (1/np.sqrt(2))* np.array([a2 - a3, a2 + a3, - a0 - a1, a0 - a1] ) ,
    (1/np.sqrt(2))*np.array([a1 - a3, - a0 + a2, - a1 - a3, a0 + a2] )
    
    ])
    
    cos_half_theta_max = 0.0
    for i in range(len(cubic_sym)):
        cos_half_theta = np.abs(cubic_sym[i][0])
        if cos_half_theta > cos_half_theta_max:
            cos_half_theta_max= cos_half_theta
        
    
    # cos_half_theta_max = 0.0
    # for i in range(len(cubic_sym)):
    #     cos_half_theta = np.abs(q_mult(cubic_sym[i], q1).dot(q2))
    #     if cos_half_theta > cos_half_theta_max:
    #                 cos_half_theta_max = cos_half_theta
    if cos_half_theta_max > 1.0:
                mis = 0.0
    else:
        mis = np.degrees( np.abs( 2 * np.arccos(cos_half_theta_max) ))
        
        
        
    return(mis)    
    

def Grain_total_curvature_surface_area_volume(mesh):
    
    # mesh_pv=pv.wrap(mesh)
    # meshfix = mf.MeshFix(mesh_pv)
    # meshfix.repair(verbose=False)
    # mesh_pv_repaired = meshfix.mesh
    # mesh=trimesh.Trimesh(mesh_pv_repaired.points,faces=mesh_pv_repaired.faces.reshape(-1, 4)[:, 1:])
    
    mesh.process()
    trimesh.repair.fill_holes(mesh)
    trimesh.repair.fix_normals(mesh, multibody=False)
    trimesh.repair.fix_inversion(mesh, multibody=False)
    # trimesh.smoothing.filter_taubin(mesh, lamb=0.5, nu=0.5, iterations=0, laplacian_operator=None)
    
    # total_mean_curvature=mesh.integral_mean_curvature/np.pi
    
    mesh_vertices=mesh.vertices
    face_adjacency_angles=mesh.face_adjacency_angles
    face_adjacency_edges=mesh.face_adjacency_edges
    # edges_unique_length=mesh.edges_unique_length
    face_adjacency_convex=mesh.face_adjacency_convex
    mean_width=0
    for i in range(0,len(face_adjacency_edges)):
        # e_i=edges_unique_length[i]
        e_i=np.linalg.norm(mesh_vertices[face_adjacency_edges[i][0]]- mesh_vertices[face_adjacency_edges[i][1]] )
        alpha_i=face_adjacency_angles[i]
        
        if face_adjacency_convex[i]==False:
            alpha_i=-alpha_i
        mean_width=mean_width+e_i*alpha_i
    mean_width=(1/(2*np.pi))*mean_width
    total_mean_curvature=mean_width
    surface_area=mesh.area
    volume=mesh.volume
    
    return(total_mean_curvature,surface_area,volume)



def Grain_imc_sa_V(mesh):
    
    # mesh_pv=pv.wrap(mesh)
    # meshfix = mf.MeshFix(mesh_pv)
    # meshfix.repair(verbose=False)
    # mesh_pv_repaired = meshfix.mesh
    # mesh=trimesh.Trimesh(mesh_pv_repaired.points,faces=mesh_pv_repaired.faces.reshape(-1, 4)[:, 1:])
    
    mesh.process()
    trimesh.repair.fill_holes(mesh)
    trimesh.repair.fix_normals(mesh, multibody=False)
    trimesh.repair.fix_inversion(mesh, multibody=False)
    # trimesh.smoothing.filter_taubin(mesh, lamb=0.5, nu=0.5, iterations=0, laplacian_operator=None)
    
    # total_mean_curvature=mesh.integral_mean_curvature/np.pi
    
    mesh_vertices=mesh.vertices
    face_adjacency_angles=mesh.face_adjacency_angles
    face_adjacency_edges=mesh.face_adjacency_edges
    # edges_unique_length=mesh.edges_unique_length
    face_adjacency_convex=mesh.face_adjacency_convex
    mean_width=0
    for i in range(0,len(face_adjacency_edges)):
        # e_i=edges_unique_length[i]
        e_i=np.linalg.norm(mesh_vertices[face_adjacency_edges[i][0]]- mesh_vertices[face_adjacency_edges[i][1]] )
        alpha_i=face_adjacency_angles[i]
        
        if face_adjacency_convex[i]==False:
            alpha_i=-alpha_i
        mean_width=mean_width+e_i*alpha_i
    mean_width=(1/(2*np.pi))*mean_width
    total_mean_curvature=mean_width
    surface_area=mesh.area
    volume=mesh.volume
    computed_properties={'imc':total_mean_curvature,'sa':surface_area,'V':volume}
    
    return(computed_properties)

def GrainBoundary_imc_sa_n_V(mesh):
    mesh_pv_repaired=pv.wrap(mesh)
    
    mesh_vertices=mesh.vertices
    face_adjacency_angles=mesh.face_adjacency_angles
    face_adjacency_edges=mesh.face_adjacency_edges
    face_adjacency_convex=mesh.face_adjacency_convex
    mean_width=0
    for i in range(0,len(face_adjacency_edges)):
        e_i=np.linalg.norm(mesh_vertices[face_adjacency_edges[i][0]]- mesh_vertices[face_adjacency_edges[i][1]] ) #edges_unique_length[i]
        alpha_i=face_adjacency_angles[i]
        
        if face_adjacency_convex[i]==False:
            alpha_i=-alpha_i
        mean_width=mean_width+e_i*alpha_i
    mean_width=(1/(2*np.pi))*mean_width
    total_mean_curvature=mean_width
    surface_area=mesh.area
    # if surface_area >0 and len(mesh.triangles)>=1:
    #     mean_curvature=np.mean(np.abs(mesh_pv_repaired.curvature()))
    #     total_mean_curvature=mesh_pv_repaired.area*mean_curvature # comment if total mean curvature from trimesh is needed
    volume=mesh.volume
    tot_curve_face=np.abs(mesh.integral_mean_curvature)
    GB_normal=np.mean(mesh.vertex_normals,axis=0).tolist()
    computed_properties={'imc':total_mean_curvature,'sa':surface_area,'GB_normal':GB_normal,'V':volume, 'imc_trimesh_positive':tot_curve_face}
    return(computed_properties)
        
def GrainBoundary_total_curvature_surface_area_volume(mesh):
    # mesh_pv=pv.wrap(mesh)
    # meshfix = mf.MeshFix(mesh_pv)
    # holes = meshfix.extract_holes()
    
    # mfix = PyTMesh(False)
    # mfix.load_array(meshfix.v, meshfix.f)
    # mfix.fill_small_boundaries(nbe=100, refine=True)
    # vert, faces = mfix.return_arrays()
    # triangles = np.empty((faces.shape[0], 4), dtype=faces.dtype)
    # triangles[:, -3:] = faces
    # triangles[:, 0] = 3

    # mesh_pv_repaired = pv.PolyData(vert, triangles)
    
    # mesh=trimesh.Trimesh(mesh_pv_repaired.points,faces=mesh_pv_repaired.faces.reshape(-1, 4)[:, 1:])
    
    # mesh.process()
    # trimesh.repair.fill_holes(mesh)
    # trimesh.repair.fix_normals(mesh, multibody=False)
    # trimesh.repair.fix_inversion(mesh, multibody=False)
    # trimesh.smoothing.filter_taubin(mesh, lamb=0.5, nu=0.5, iterations=0, laplacian_operator=None)
    
    # total_mean_curvature=mesh.integral_mean_curvature/np.pi
    mesh_pv_repaired=pv.wrap(mesh)
    
    mesh_vertices=mesh.vertices
    face_adjacency_angles=mesh.face_adjacency_angles
    face_adjacency_edges=mesh.face_adjacency_edges
    face_adjacency_convex=mesh.face_adjacency_convex
    mean_width=0
    for i in range(0,len(face_adjacency_edges)):
        e_i=np.linalg.norm(mesh_vertices[face_adjacency_edges[i][0]]- mesh_vertices[face_adjacency_edges[i][1]] ) #edges_unique_length[i]
        alpha_i=face_adjacency_angles[i]
        
        if face_adjacency_convex[i]==False:
            alpha_i=-alpha_i
        mean_width=mean_width+e_i*alpha_i
    mean_width=(1/(2*np.pi))*mean_width
    total_mean_curvature=mean_width
    surface_area=mesh.area
    # if surface_area >0 and len(mesh.triangles)>=1:
    #     mean_curvature=np.mean(np.abs(mesh_pv_repaired.curvature()))
    #     total_mean_curvature=mesh_pv_repaired.area*mean_curvature # comment if total mean curvature from trimesh is needed
    volume=mesh.volume
    tot_curve_face=np.abs(mesh.integral_mean_curvature)
    
    return(total_mean_curvature,surface_area,volume,tot_curve_face)

     

def GrainBoundary_total_curvature_surface_area_volume_extended(mesh_list):
    
    m0=mesh_list[0]
    m1=mesh_list[1]
    m0_total_mean_curvature,m0_surface_area,m0_volume=GrainBoundary_total_curvature_surface_area_volume(m0)
    m1_total_mean_curvature,m1_surface_area,m1_volume=GrainBoundary_total_curvature_surface_area_volume(m1)

    if (m0_total_mean_curvature!=0 and m0_surface_area!=0 and m1_total_mean_curvature!=0 and m1_surface_area!=0 ):
        surface_area=0.5*(np.abs(m0_surface_area)+np.abs(m1_surface_area))
        
        total_mean_curvature=0.5*(np.abs(m0_total_mean_curvature)+np.abs(m1_total_mean_curvature))
        
        mean_curvature=0.5*(np.abs(m0_total_mean_curvature/m0_surface_area) + np.abs(m1_total_mean_curvature/m1_surface_area))
        
    elif (m0_total_mean_curvature==0 or  m0_surface_area==0 or m1_total_mean_curvature==0 or m1_surface_area==0 ):
        total_mean_curvature=np.nan
        surface_area=np.nan
        mean_curvature=np.nan
    
    return(total_mean_curvature,mean_curvature,surface_area)
    
    
        
        
        
        
    






def Microstructure_topology(InnerAtom,GrainBoundary,TripleLine,grain_property_current,grain_boundary_property_current,grain_boundary_property_extended_current,triple_line_length_current,timestep,q0,q1,q2,q3,grain_mapping_data,X,GrainId_list,EdgeLength):
    
    '''
    This function aggregates relevant information for grains, grain boundaries and triple junctions for further anlysis
    '''
    GrainBoundary_misorientation=[]
    GrainBoundary_number_of_triple_line=[]
    GrainBoundary_length_of_triple_line=[]
    
    
    
    for i in range(0,len(GrainBoundary)):
    
        Grain_Boundary_id=GrainBoundary[i][0]
        q1_i,q2_i=misorientation_quaternions(timestep,Grain_Boundary_id, q0,q1,q2,q3,grain_mapping_data)
        mis_i=cubically_equivalent_quaternions (q1_i,q2_i)
        
        length_of_triple_line_associated_with_grain_boundary=[j[1] for j in triple_line_length_current if (Grain_Boundary_id[0] in j[0] and Grain_Boundary_id[1] in j[0])]
        total_length_of_triple_line_associated_with_grain=sum(length_of_triple_line_associated_with_grain_boundary)
        
        
        triple_lines_associated_with_grain_boundary=[j[0] for j in  TripleLine if  (Grain_Boundary_id[0] in j[0] and Grain_Boundary_id[1] in j[0]) ]
        number_of_triple_lines_associated_with_grain_boundary=len(triple_lines_associated_with_grain_boundary)
        
        GrainBoundary_number_of_triple_line.append([Grain_Boundary_id,number_of_triple_lines_associated_with_grain_boundary, triple_lines_associated_with_grain_boundary])
        GrainBoundary_misorientation.append([Grain_Boundary_id,mis_i,triple_lines_associated_with_grain_boundary,[q1_i,q2_i]])
        GrainBoundary_length_of_triple_line.append([Grain_Boundary_id,total_length_of_triple_line_associated_with_grain,triple_lines_associated_with_grain_boundary])
        
        
    GrainBoundary_total_curvature_mean_curvature_surface_area_current=[]  
    ##############################################
    # for i in range (0,len(grain_boundary_property_current)):
    #     Grain_Boundary_id=grain_boundary_property_current[i][0]
    #     mesh_list=grain_boundary_property_current[i][1]
    #     # total_mean_curvature,surface_area,volume=GrainBoundary_total_curvature_surface_area_volume(mesh)
    #     total_mean_curvature,mean_curvature,surface_area=GrainBoundary_total_curvature_surface_area_volume_extended(mesh_list)
    #     GrainBoundary_total_curvature_mean_curvature_surface_area_current.append([Grain_Boundary_id,total_mean_curvature,mean_curvature,surface_area])
        
    ############################################################
    for i in range (0,len(grain_boundary_property_extended_current)):
        Grain_Boundary_id=grain_boundary_property_extended_current[i][0]
        mesh=grain_boundary_property_extended_current[i][1]
        total_mean_curvature,surface_area,volume,tot_curve_face=GrainBoundary_total_curvature_surface_area_volume(mesh)
        if surface_area !=0:
            mean_curvature=total_mean_curvature/abs(surface_area)
            
        else:
            mean_curvature=np.nan
        GrainBoundary_total_curvature_mean_curvature_surface_area_current.append([Grain_Boundary_id,total_mean_curvature,mean_curvature,surface_area,tot_curve_face])
        
    
    #############################################################    
    Grain_number_of_grain_boundaries=[]
    Grain_average_misorientation=[]
    
    for i in range(0,len(InnerAtom)):
    
        Grain_id=InnerAtom[i][0]
        misorientation_of_gb_associated_with_grain=[j[1] for j in GrainBoundary_misorientation if (Grain_id[0] in j[0])]
        
        # area_of_gb_associated_with_grain=[j[3] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_current if (Grain_id[0] in j[0])]
        
        # wighted_area_of_gb_associated_with_grain=np.array(area_of_gb_associated_with_grain)/np.nansum(np.array(area_of_gb_associated_with_grain))
        
        # area_wighted_average_misorientation_of_gb_associated_with_grain=np.array(misorientation_of_gb_associated_with_grain)*wighted_area_of_gb_associated_with_grain
        average_misorientation_of_gb_associated_with_grain= sum(misorientation_of_gb_associated_with_grain)/len(misorientation_of_gb_associated_with_grain) #np.nanmean(area_wighted_average_misorientation_of_gb_associated_with_grain)
        grain_boundary_associated_with_grain=[j[0] for j in  GrainBoundary if  (Grain_id[0] in j[0]) ]
        
        number_of_grain_boundaries_associated_with_grain=len(grain_boundary_associated_with_grain)
        
        
        Grain_number_of_grain_boundaries.append([Grain_id,number_of_grain_boundaries_associated_with_grain, grain_boundary_associated_with_grain])
        Grain_average_misorientation.append([Grain_id,average_misorientation_of_gb_associated_with_grain,grain_boundary_associated_with_grain] )
        
        
    Grain_number_of_triple_lines=[]
    Grain_length_of_triple_lines=[]
    
    for i in range(0,len(InnerAtom)):
    
        Grain_id=InnerAtom[i][0]
        length_of_triple_line_associated_with_grain=[j[1] for j in triple_line_length_current if (Grain_id[0] in j[0])]
        lam_i_of_triple_line_associated_with_grain=[j[3] for j in triple_line_length_current if (Grain_id[0] in j[0])]
        average_of_lam_i_associated_with_grain=np.nanmean(lam_i_of_triple_line_associated_with_grain)
        total_length_of_triple_line_associated_with_grain=sum(length_of_triple_line_associated_with_grain)
        triple_line_associated_with_grain=[j[0] for j in  TripleLine if  (Grain_id[0] in j[0]) ]
        number_of_triple_lines_associated_with_grain=len(triple_line_associated_with_grain)
        
        Grain_number_of_triple_lines.append([Grain_id,number_of_triple_lines_associated_with_grain, triple_line_associated_with_grain])
        Grain_length_of_triple_lines.append([Grain_id,total_length_of_triple_line_associated_with_grain, triple_line_associated_with_grain,average_of_lam_i_associated_with_grain])
   
    grain_total_curvature_surface_area_volume_current=[]
    grain_total_curvature_surface_area_volume_from_face_current=[]
    
    for i in range (0,len(grain_property_current)):
        Grain_id=grain_property_current[i][0]
        mesh=copy.deepcopy(grain_property_current[i][1]) #Grain_mesh_ovito(GrainId_list,X,EdgeLength,Grain_id)
        total_mean_curvature,surface_area,volume=Grain_total_curvature_surface_area_volume(mesh)
        grain_total_curvature_surface_area_volume_current.append([Grain_id,total_mean_curvature,surface_area,volume])
        
        grain_boundary_mesh_associated_with_grain=[j[1] for j in  grain_boundary_property_extended_current if  (Grain_id == j[0][0]) ]
        tsv=[GrainBoundary_total_curvature_surface_area_volume(copy.deepcopy(j)) for j in  grain_boundary_mesh_associated_with_grain]
        
        total_mean_curvature_from_face= sum([j[0] for j in  tsv])
        
        surface_area_from_face= sum([j[1] for j in  tsv])
        
        volume_from_face= sum([j[2] for j in  tsv])
        
        total_mean_curvature_from_face_positive= sum([j[3] for j in  tsv])
        grain_total_curvature_surface_area_volume_from_face_current.append([Grain_id,total_mean_curvature_from_face,surface_area_from_face,volume_from_face,total_mean_curvature_from_face_positive])
    
    return(Grain_number_of_grain_boundaries, Grain_number_of_triple_lines,Grain_average_misorientation,Grain_length_of_triple_lines,grain_total_curvature_surface_area_volume_current, GrainBoundary_number_of_triple_line,GrainBoundary_misorientation,GrainBoundary_length_of_triple_line,GrainBoundary_total_curvature_mean_curvature_surface_area_current,grain_total_curvature_surface_area_volume_from_face_current)
      