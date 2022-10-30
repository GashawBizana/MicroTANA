import trimesh
import pyvista as pv
import numpy as np
import time
import pymeshfix as mf
import matplotlib as mpl
import matplotlib.pylab as plt
import copy


def distance_from_von_neuman(grain_id_k_i_initial,grain_boundary_k_i,grain_boundary_k_i_triangleId,Grain_VonNeuman_Data_step_k):

    
    dis=[j[14] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[0]== j[0] ][0]
    vol1=[j[7] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[0]== j[0] ][0]
    
    vol2=[j[7] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[1]== j[0] ]
    
    G_tcf1=[j[15] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[0]== j[0] ][0]
    
    G_tcf2=[j[15] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[1]== j[0] ]
    
    if not vol2:
        
        vol2=0.0
        G_tcf2=0.0
    else:
        vol2=vol2[0]
        G_tcf2=G_tcf2[0]
    # vol1=min(vol1,vol2)
    # vol2=min(vol1,vol2)
    triangle_distance=[]
    for i in range(0,len(grain_boundary_k_i_triangleId)):
        
       
        
        dis_tri=1/9*(dis[grain_boundary_k_i_triangleId[i][0]] + dis[grain_boundary_k_i_triangleId[i][1]]+ dis[grain_boundary_k_i_triangleId[i][2]])
        
        triangle_distance.append(dis_tri)
    distance=np.nanmean(triangle_distance)
    
    curv=np.abs(grain_boundary_k_i.integral_mean_curvature/grain_boundary_k_i.area)
    
    return(distance,curv,vol1,vol2,G_tcf1,G_tcf2)
        

def distanace_between_current_and_next(grain_boundary_step_i, grain_boundary_step_i_plus_1,EdgeLength_step_i_plus_1,Grain_VonNeuman_Data_step_k,grain_id_k_i_initial,grainMesh):
    
    vol1=[j[7] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[0]== j[0] ][0]
    
    vol2=[j[7] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[1]== j[0] ]
    
    G_tcf1=[j[15] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[0]== j[0] ][0]
    
    G_tcf2=[j[15] for j in Grain_VonNeuman_Data_step_k if grain_id_k_i_initial[1]== j[0] ]
    if not vol2:
        
        vol2=0.0
        G_tcf2=0.0
    else:
        vol2=vol2[0]
        G_tcf2=G_tcf2[0]
    
    
    grain_boundary_step_i_plus_1=GrainsRecenter(grain_boundary_step_i,grain_boundary_step_i_plus_1,EdgeLength_step_i_plus_1)
    
    grain_boundary_step_i_plus_1.remove_unreferenced_vertices()
    grain_boundary_step_i.remove_unreferenced_vertices()
    # cloud_grain_boundary_step_i= pv.PolyData(grain_boundary_step_i)
    # cloud_grain_boundary_step_i_plus_1=pv.PolyData(grain_boundary_step_i_plus_1)
    
    
    
    # h0= cloud_grain_boundary_step_i.delaunay_2d()
    # h1 = cloud_grain_boundary_step_i_plus_1.delaunay_2d()
    
    h0=pv.wrap(grain_boundary_step_i)
    h0=h0.subdivide_adaptive(max_n_passes=1)
    h1=pv.wrap(grain_boundary_step_i_plus_1)
    h1=h1.subdivide_adaptive(max_n_passes=1)
    # distance=np.sqrt(np.sum((np.array(h0.center) - np.array(h1.center)) ** 2))
    def array_dis(h0,h1):
        h0n = h0.compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=True)
        
        h0n["distances"] = np.empty(h0.n_points)
        for i in range(h0n.n_points):
            p = h0n.points[i]
            vec = h0n["Normals"][i] * h0n.length
            p0 = p - vec
            p1 = p + vec
            ip, ic = h1.ray_trace(p0, p1, first_point=True)
            dist = np.sqrt(np.sum((ip - p) ** 2))
            h0n["distances"][i] = dist
        
        # Replace zeros with nans
        mask = h0n["distances"] == 0
        h0n["distances"][mask] = np.nan
        hh=h0n["distances"]
        distance=np.mean(hh[hh<=7])
        # distance=np.nanmean( np.where(h0n["distances"] <= 7, np.nan, h0n["distances"]))
        return(distance)
    
    def exact_dis(h0,h1):
        closest_cells, closest_points = h1.find_closest_cell(h0.points, return_closest_point=True)
        d_exact = np.linalg.norm(h0.points - closest_points, axis=1)
        h0["distances"] = d_exact
        distance=np.mean(np.abs(d_exact))
        
        return(distance)
    
    def ray_dis(m1,m2):
        
        t=trimesh.proximity.longest_ray(m2, m1.vertices, -1*m1.vertex_normals)
        t1=trimesh.proximity.longest_ray(m2, m1.vertices, m1.vertex_normals)
        tt=[]
        for i in range(0,len(t)):
            dt=t[i]
            dt1=t1[i]
            if dt==np.inf and dt1==np.inf:
                tt.append(np.nan)
            elif dt <=dt1:
                tt.append(-dt)
            else:
                tt.append(dt1)
        
        dis=np.nanmean(tt)
        return(dis)
    
    def volume_dis(h0,h1):
        
        distance=np.abs(np.abs(h0.volume)-np.abs(h1.volume))/(0.5*np.abs(np.abs(h0.area)+np.abs(h1.area)))
        
        return(distance)
    
    def centroid_dis(h0,h1):
        
        distance=np.linalg.norm(np.array(h0.center_of_mass())-np.array(h1.center_of_mass()))
        return(distance)
    
    # dis1=array_dis(h0,h1)
    # dis2=array_dis(h1,h0)
    
    dis1=exact_dis(h0,h1)
    dis2=exact_dis(h1,h0)
    
    # dis1=volume_dis(h0,h1)
    # dis2=volume_dis(h1,h0)
    
    # dis1=ray_dis(grain_boundary_step_i,grain_boundary_step_i_plus_1)
    dis1=centroid_dis(h1,h0)
    
    distance=0.5*(dis1+dis1)
    curv=np.mean(h0.curvature())
                   
    
    # mean=np.nanmean(h0n["distances"])
    # std=np.nanstd(h0n["distances"])
    # distance= np.nanmean(h0n["distances"][h0n["distances"]  <= 100])
    # distance= np.nanmean(h0n["distances"])
    # distance=np.nanmean( np.where(h0n["distances"] <= 2, np.nan, h0n["distances"]))
    
    # closest_cells, closest_points = h1.find_closest_cell(h0.points, return_closest_point=True)
    # d_exact = np.linalg.norm(h0.points - closest_points, axis=1)
    # h0["distances"] = d_exact
    # distance=np.mean(d_exact)
    
    
    return(distance,curv,vol1,vol2,G_tcf1,G_tcf2)

def Initial_Grain_Id(grain_mapping_data, grain_boundary_forming_grains_n,time_step_n):
    n0=grain_boundary_forming_grains_n[0]
    n1=grain_boundary_forming_grains_n[1]
    grain_mapping_data_time_step_n=grain_mapping_data[time_step_n]
    index_n0=np.where(grain_mapping_data_time_step_n==n0)
    index_n1=np.where(grain_mapping_data_time_step_n==n1)
    
    return([int(index_n0[0]),int(index_n1[0])])

def GrainIdMapper(grain_mapping_data, grain_boundary_forming_grains_n,time_step_n,step):
    
    n0=grain_boundary_forming_grains_n[0]
    n1=grain_boundary_forming_grains_n[1]
    grain_mapping_data_time_step_n=grain_mapping_data[time_step_n]
    index_n0=np.where(grain_mapping_data_time_step_n==n0)
    index_n1=np.where(grain_mapping_data_time_step_n==n1)
    grain_mapping_data_time_step_m=grain_mapping_data[time_step_n+step]
    
    if len(index_n0[0])==0:
        m0=np.nan
    else:
        m0=(grain_mapping_data_time_step_m[index_n0])[0]
    if len(index_n1[0])==0:
        m1=np.nan
    else:
        m1=(grain_mapping_data_time_step_m[index_n1])[0]
    
    if np.isnan(m0) or np.isnan(m1):
        
        grain_boundary_forming_grains_m=[m0,m1]
    else:
        
        grain_boundary_forming_grains_m=[int(m0),int(m1)]
    
    return(grain_boundary_forming_grains_m)

def GrainRepair(mesh):
    
    mesh_pv=pv.wrap(mesh)
    meshfix = mf.MeshFix(mesh_pv)
    meshfix.repair(verbose=False)
    mesh_pv_repaired = meshfix.mesh
    
    mesh=trimesh.Trimesh(mesh_pv_repaired.points,faces=mesh_pv_repaired.faces.reshape(-1, 4)[:, 1:])
    
    mesh.process(validate=True, merge_tex=None, merge_norm=None)
    trimesh.repair.fill_holes(mesh)
    trimesh.repair.fix_normals(mesh, multibody=False)
    trimesh.repair.fix_inversion(mesh, multibody=False)
    
    return(mesh)

def GrainsRecenter(grain_mesh_n_0,grain_mesh_m_1,EdgeLength_step_i_plus_1):
    
    grain_mesh_n_0_cm=np.mean(grain_mesh_n_0.vertices,axis=0)
    grain_mesh_m_1_cm=np.mean(grain_mesh_m_1.vertices,axis=0)
    
    diff_center=np.rint((grain_mesh_n_0_cm-grain_mesh_m_1_cm)/EdgeLength_step_i_plus_1)    
    periodic_translation_factor=diff_center*EdgeLength_step_i_plus_1
    grain_mesh_m_1.vertices= grain_mesh_m_1.vertices + periodic_translation_factor
    
    return(grain_mesh_m_1)

def IntersectionVolume(grain_mesh_n_0, grain_mesh_m_1,EdgeLength_step_k_plus_1):
    
    grain_mesh_m_1=GrainsRecenter(grain_mesh_n_0, grain_mesh_m_1,EdgeLength_step_k_plus_1)
    
    # grain_mesh_n_0=GrainRepair(grain_mesh_n_0)
    
    # grain_mesh_m_1=GrainRepair(grain_mesh_m_1)
    
    Intersection=trimesh.boolean.intersection([grain_mesh_n_0,grain_mesh_m_1], engine='blender')
    # Intersection = grain_mesh_n_0.boolean_intersection(grain_mesh_m_1)
    
    return(Intersection)



def GrainBoundary_distance_computation(EdgeLength_for_all_timesteps,grain_boundary_property_for_all_timesteps,GrainBoundary_length_of_triple_line_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps,GrainBoundary_number_of_triple_line_for_all_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps,grain_property_for_all_timesteps,grain_boundary_property_extended_for_all_timesteps,grain_mapping_data,step,Grain_VonNeuman_Data_for_timesteps):
    
    Grain_Boundary_Velocity_Data_for_all_timesteps=[]
    #len(grain_boundary_property_for_all_timesteps)-step
    #len(grain_boundary_property_for_all_timesteps)-step
    for k in range(0,len(grain_boundary_property_for_all_timesteps)-step):
        print(k)
        Grain_VonNeuman_Data_step_k=Grain_VonNeuman_Data_for_timesteps[k]
        
        
        EdgeLength_step_k= EdgeLength_for_all_timesteps[k]
        EdgeLength_step_k_plus_1=EdgeLength_for_all_timesteps[k+step]
        # print(EdgeLength_step_k_plus_1)
        
        grain_boundary_property_step_k=grain_boundary_property_extended_for_all_timesteps[k]
        grain_boundary_property_step_k_plus_1=grain_boundary_property_extended_for_all_timesteps[k+step]
        
        grain_boundary_total_curvature_step_k=[j[1] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]]
        grain_boundary_total_surface_area_step_k=[j[3] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]]
        grain_boundary_mean_curvature_step_k =[j[2] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]]
        
        
        # grain_boundary_forming_step_k=[j[0] for j in grain_boundary_property_step_k]
        
        grain_boundary_forming_step_k=[j[0] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]]
        
        grain_property_step_k=grain_property_for_all_timesteps[k]
        grain_property_step_k_plus_1=grain_property_for_all_timesteps[k+step]
        
        GrainBoundary_length_of_triple_line_step_k=GrainBoundary_length_of_triple_line_for_all_timesteps[k]
        GrainBoundary_length_of_triple_line_step_k_plus_1=GrainBoundary_length_of_triple_line_for_all_timesteps[k+step]
        GrainBoundary_misorientation_step_k=GrainBoundary_misorientation_for_all_timesteps[k]
        GrainBoundary_misorientation_step_k_plus_1=GrainBoundary_misorientation_for_all_timesteps[k+step]
        GrainBoundary_number_of_triple_line_step_k=GrainBoundary_number_of_triple_line_for_all_timesteps[k]
        GrainBoundary_number_of_triple_line_step_k_plus_1=GrainBoundary_number_of_triple_line_for_all_timesteps[k+step]
        
        Grain_Boundary_Velocity_Data_k=[]
        IntersectionData=[]
        for i in range(0,len(grain_boundary_total_curvature_step_k)):
            
            #print(i)
            grain_boundary_total_curvature_step_k_i=grain_boundary_total_curvature_step_k[i]
            grain_boundary_total_surface_area_step_k_i=grain_boundary_total_surface_area_step_k[i]
            grain_boundary_mean_curvature_step_k_i =grain_boundary_mean_curvature_step_k[i]
            
            # if grain_boundary_total_surface_area_step_k_i==0:
            #     grain_boundary_mean_curvature_step_k_i=0
            # elif grain_boundary_total_surface_area_step_k_i !=0:
            #     grain_boundary_mean_curvature_step_k_i=grain_boundary_total_curvature_step_k_i/grain_boundary_total_surface_area_step_k_i
            
            
            grain_id_k_i=grain_boundary_forming_step_k[i]
            #print(grain_id_k_i)
            grain_id_k_i_initial=Initial_Grain_Id(grain_mapping_data, grain_id_k_i ,k)
            
            grain_id_k_plus_1_i=GrainIdMapper(grain_mapping_data, grain_id_k_i ,k,step)
            
            # print(grain_id_k_plus_1_i)
            grain_boundary_k_i_triangleId=[j[2] for j in grain_boundary_property_step_k if (grain_id_k_i[0] == j[0][0] and grain_id_k_i[1] == j[0][1])]
            
            grain_boundary_k_i=[j[1] for j in grain_boundary_property_step_k if (grain_id_k_i[0] == j[0][0] and grain_id_k_i[1] == j[0][1])]
            
            grain_k_i=[j[1] for j in grain_property_step_k if j[0]== grain_id_k_i[0]]
            grain_k_i_id=[j[0] for j in grain_property_step_k if j[0]== grain_id_k_i[0]]
            # print(grain_k_i_id)
            
            grain_boundary_k_plus_1_i=[j[1] for j in grain_boundary_property_step_k_plus_1 if (grain_id_k_plus_1_i[0] == j[0][0] and grain_id_k_plus_1_i[1] == j[0][1])]
            
            
            grain_id_k_i0=grain_id_k_i[0]
            grain_id_k_i1=grain_id_k_i[1]
            grain_id_k_plus_1_i0=grain_id_k_plus_1_i[0]
            grain_id_k_plus_1_i1=grain_id_k_plus_1_i[1]
            
            
            GrainBoundary_length_of_triple_line_step_k_i=[j[1] for j in GrainBoundary_length_of_triple_line_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            GrainBoundary_misorientation_step_k_i= [j[1] for j in GrainBoundary_misorientation_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            GrainBoundary_quat_step_k_i= [j[3] for j in GrainBoundary_misorientation_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            GrainBoundary_number_of_triple_line_step_k_i= [j[1] for j in GrainBoundary_number_of_triple_line_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            
            if (np.isnan(grain_id_k_plus_1_i).any()) ==False and ((not grain_boundary_k_plus_1_i)== False) and (len(grain_boundary_k_i[0].faces) >=1) and len(grain_boundary_k_plus_1_i[0].faces) >=1:

                print(grain_id_k_i)
                
                grain_mesh_k_i0L=[j[1] for j in grain_property_step_k if j[0]==grain_id_k_i0]
                grain_mesh_k_i1L=[j[1] for j in grain_property_step_k if j[0]==grain_id_k_i1]
                grain_mesh_k_plus_1_i1L=[j[1] for j in grain_property_step_k_plus_1 if j[0]==grain_id_k_plus_1_i1]
                grain_mesh_k_plus_1_i0L=[j[1] for j in grain_property_step_k_plus_1 if j[0]==grain_id_k_plus_1_i0]
                
                if (not grain_mesh_k_i0L) or (not grain_mesh_k_i1L) or (not grain_mesh_k_plus_1_i1L) or (not grain_mesh_k_plus_1_i0L):
                    continue
                
                grain_mesh_k_i0=grain_mesh_k_i0L[0]
                grain_mesh_k_i1=grain_mesh_k_i1L[0]
                grain_mesh_k_plus_1_i1=grain_mesh_k_plus_1_i1L[0]
                grain_mesh_k_plus_1_i0=grain_mesh_k_plus_1_i0L[0]
                
                Extract_computed_intersection=[val[1] for val in IntersectionData if (grain_id_k_i1==val[0][0] and grain_id_k_i0==val[0][1])]
                #print(Extract_computed_intersection)
                if (not Extract_computed_intersection)==True:
                    #print('No')
                
                    Intersection0=IntersectionVolume(copy.deepcopy(grain_mesh_k_i0), copy.deepcopy(grain_mesh_k_plus_1_i1),EdgeLength_step_k_plus_1)
                    Intersection1=IntersectionVolume(copy.deepcopy(grain_mesh_k_i1), copy.deepcopy(grain_mesh_k_plus_1_i0),EdgeLength_step_k_plus_1)
                    
                    if ((str(type(Intersection0)) == "<class 'trimesh.scene.scene.Scene'>")):
                        IntVol_0=0
                    else:
                        IntVol_0=Intersection0.volume
                        
                    if ((str(type(Intersection1)) == "<class 'trimesh.scene.scene.Scene'>")):
                        IntVol_1=0
                    else:
                        IntVol_1=Intersection1.volume
                    
                    
                    average_area= (Intersection0.area + Intersection1.area)/2
                        # if Intersection0.area !=0:
                elif (not Extract_computed_intersection)==False:
                    #print('yes')
                    IntVol_0=Extract_computed_intersection[0][1]
                    IntVol_1=Extract_computed_intersection[0][0]
                    average_area=Extract_computed_intersection[0][2]
                        
                    
                IntersectionData.append([grain_id_k_i,[IntVol_0,IntVol_1,average_area]])
                # velocity,curv_py,vol1,vol2=distanace_between_current_and_next(copy.deepcopy(grain_boundary_k_i[0]),copy.deepcopy(grain_boundary_k_plus_1_i[0]),EdgeLength_step_k_plus_1,Grain_VonNeuman_Data_step_k,grain_id_k_i_initial,grain_k_i[0])
                # print(velocity)
                velocity,curv_py,vol1,vol2,G_tcf1,G_tcf2=distance_from_von_neuman(grain_id_k_i_initial,grain_boundary_k_i[0],grain_boundary_k_i_triangleId[0],Grain_VonNeuman_Data_step_k)
                
                    # velocity=2*np.abs(Intersection0.volume)/(abs(grain_boundary_k_i[0].area)+abs(grain_boundary_k_plus_1_i[0].area))
                Grain_Boundary_Velocity_Data_k.append([grain_id_k_i,grain_id_k_i_initial,grain_id_k_plus_1_i,velocity,grain_boundary_total_curvature_step_k_i,grain_boundary_mean_curvature_step_k_i,grain_boundary_total_surface_area_step_k_i,np.nan,GrainBoundary_length_of_triple_line_step_k_i,GrainBoundary_number_of_triple_line_step_k_i,GrainBoundary_misorientation_step_k_i,curv_py,vol1,vol2,grain_boundary_k_i[0],G_tcf1,G_tcf2,IntVol_0,IntVol_1,average_area,GrainBoundary_quat_step_k_i])
        
        '''
        0- grain_id_k_i, 1- grain_id_k_i_initial, 2- grain_id_k_plus_1_i, 3- Velocity, 4- grain_boundary_total_curvature_step_k_i, 5-grain_boundary_mean_curvature_step_k_i, 6- grain_boundary_total_surface_area_step_k_i, 7- grain_boundary_total_volume_step_k_i, 8- GrainBoundary_length_of_triple_line_step_k_i, 9- GrainBoundary_number_of_triple_line_step_k_i, 10- GrainBoundary_misorientation_step_k_i
        '''
        Grain_Boundary_Velocity_Data_for_all_timesteps.append(Grain_Boundary_Velocity_Data_k)
    
    return(Grain_Boundary_Velocity_Data_for_all_timesteps)

step=1
t = time.time()
Grain_Boundary_Velocity_Data_for_all_timesteps=GrainBoundary_distance_computation(EdgeLength_for_all_timesteps,grain_boundary_property_for_all_timesteps,GrainBoundary_length_of_triple_line_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps,GrainBoundary_number_of_triple_line_for_all_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps,grain_property_for_all_timesteps,grain_boundary_property_extended_for_all_timesteps,grain_mapping_data,step,Grain_VonNeuman_Data_for_timesteps)


curvlimit=[-10000,100000]
dislimit=[-100000,100000]
mislimit=[0,70]
arealimit=[]
grain1=[0,150]
grain2=[0,150]

# Grain_Boundary_Velocity_Data_for_all_timesteps=Grain_Boundary_Velocity_Data_for_all_timesteps_ray

dis,curv,mis,tot_cur,area,GB_ltl,GB_ntl,curv_py,GB_v1,GB_v2,gId,ts,GB_tcf1,GB_tcf2,Int_vol1,Int_vol2,Int_a,gbq,GB_n=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
LE=len(Grain_Boundary_Velocity_Data_for_all_timesteps)
for i in range(0,LE):
    dis0=[k[3] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1]) ]
    
    curv0=[k[5] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    mis0=[k[10] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    GB_ltl0=[k[8] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    GB_ntl0=[k[9] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    tot_cur0=[k[4] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    area0=[k[6] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    curv_py0=[k[11] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    GB_v10=[k[12] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    GB_v20=[k[13] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    GB_n0=[np.mean(k[14].vertex_normals,axis=0).tolist()for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    GB_tcf10=[k[15] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    GB_tcf20=[k[16] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    Int_vol10=[k[17] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    Int_vol20=[k[18] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    Int_a0=[k[19] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    gbq0=[np.ravel(k[20]).tolist() for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    
    gId0=[k[1] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]
    ts0=[i for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[0]).any()) ==False and dislimit[0] < k[3] < dislimit[1] and curvlimit[0] < k[5] <curvlimit[1] 
and mislimit[0] < k[10] < mislimit[1] and (grain1[0]<=k[1][0]<=grain1[1] and grain2[0]<=k[1][1]<=grain2[1])]


    
    
    dis=dis+dis0
    curv=curv+curv0
    mis=mis+mis0
    tot_cur=tot_cur+tot_cur0
    area=area+area0
    GB_ntl=GB_ntl+GB_ntl0
    GB_ltl=GB_ltl+GB_ltl0
    curv_py=curv_py0+curv_py
    GB_v1=GB_v1+GB_v10
    GB_v2=GB_v2+GB_v20
    GB_n=GB_n+GB_n0
    GB_tcf1=GB_tcf1+GB_tcf10
    GB_tcf2=GB_tcf2+GB_tcf20
    Int_vol1=Int_vol1+Int_vol10
    Int_vol2=Int_vol2+Int_vol20
    Int_a=Int_a+Int_a0
    gbq=gbq+gbq0
    gId=gId+gId0
    ts=ts+ts0




dt=1/(step)
view=np.array(mis)
# xa=np.array(curv_py)
xa=np.power(np.array(GB_v1)*0.001,1/3)
xa=np.abs(np.array(tot_cur))
# xa=np.power(np.array(GB_v2)*0.001,1/3)
# ya=np.power(np.array(GB_v2)*0.001,1/3)
# xa=0.5*np.power(np.array(GB_v2)*0.001,1/3)+ 0.5*np.power(np.array(GB_v1)*0.001,1/3)
# ya=xa
#xa=np.array(dis)*0.1*dt
ya=((-np.array(Int_vol1)+np.array(Int_vol2))/(Int_a))*0.1*dt
# xa=view
#xa=np.abs(-np.array(tot_cur)*0.1/(np.array(area)*0.01))
# xa=np.abs(np.array(area)*0.1)
# xa=np.array(GB_v2)*0.001+ np.array(GB_v1)*0.001
# xa=np.power(xa,-1/3)
# xa=np.array(curv_py)
cmaplist=[cmap(i) for i in range(cmap.N)]
cmaplist[0]=(0.5,0.5,0.5,1.0)
cmap=mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist,cmap.N)

# bounds=np.linspace(np.rint(min(mis)),np.rint(max(mis)),2*int(np.rint(max(mis))-np.rint(min(mis))))
bounds=np.linspace(min(view),max(view),256)
norm=mpl.colors.BoundaryNorm(bounds, cmap.N)
# coef = np.polyfit(xa,ya,1)
# poly1d_fn = np.poly1d(coef)
# xnew=np.linspace(min(xa),max(xa),10)
fig = plt.figure()
ax = fig.add_subplot()

sc=ax.scatter(xa,ya,c=view,cmap=cmap,norm=norm,s=25,marker='+')

# ax.set_xlabel('Total curvature ($\AA$)')
ax.set_xlabel('curvature (1/nm)',fontweight='bold',size=12)
# ax.set_ylabel('phi')
ax.set_ylabel('Velocity (nm/ns)',fontweight='bold',size=12)
cb=plt.colorbar(sc)
cb.set_label('Disorientation (deg)',fontweight='bold',size=12)
# plt.plot(xnew,poly1d_fn(xnew))

plt.show()

#Grain_Boundary_Velocity_Data_for_all_timesteps[0][18][4]

dict = {'dis': dis, 'curv': curv, 'mis': mis,'tot_cur':tot_cur,'area':area,'GB_ntl':GB_ntl,'GB_ltl':GB_ltl,'curv_py':curv_py,'GB_v1':GB_v1,'GB_v2':GB_v2,'GB_tcf1':GB_tcf1,'GB_tcf2':GB_tcf2,'Int_vol1':Int_vol1,'Int_vol2':Int_vol2,'Int_a':Int_a,'gbq':gbq,'gId':gId,'ts':ts,'GB_n':GB_n}

df = pd.DataFrame(dict)
import os  
os.makedirs('E:/MicroTANA/GB_data_new1', exist_ok=True)  
df.to_csv('E:/MicroTANA/GB_data_new1/GB_all_new1.csv')

# import matplotlib.pyplot as plt
# plt.rcParams.update({'figure.figsize':(8,6), 'figure.dpi':100})
# plt.hist(np.array(mis), bins=None,density = True, color = "blue", label = "Input Microstructure",histtype='bar', ec='black')


# x=np.linspace(0,62.8,100)

# y=[]

# for i in range(0,len(x)):
#     if 0<=x[i]<=45:
#         yi=(2/15)*(1 - np.cos(np.radians(x[i])))
#     elif 45<x[i]<=60:
#         yi=(2/15)*(3*(np.sqrt(2)-1)*np.sin(np.radians(x[i]))-2*(1-np.cos(np.radians(x[i]))))
#     elif 60<x[i]<=60.72:
#         yi= (2/15)*((3*(np.sqrt(2)-1) + 4/np.sqrt(3))*np.sin(np.radians(x[i]))-6*(1-np.cos(np.radians(x[i]))))
#     elif 60.72<=x[i]<=62.8:
        
#         zzz=(2/15)*((3*(np.sqrt(2)-1) + 4/np.sqrt(3))*np.sin(np.radians(x[i]))-6*(1-np.cos(np.radians(x[i]))))
#         xx=(np.sqrt(2)-1) / (1-((np.sqrt(2)-1)**(2))*((np.tan(0.5*np.radians(x[i])))**(-2)))**(1/2)
        
#         yy=((np.sqrt(2)-1)**2) / (3-(np.tan(0.5*np.radians(x[i])))**(-2))**(1/2)
        
#         xxx=-(8/(5*np.pi))*(2*(np.sqrt(2)-1)*np.arccos(xx*(np.tan(0.5*np.radians(x[i])))**(-1)) + (1/np.sqrt(3))*np.arccos(yy*(np.tan(0.5*np.radians(x[i])))**(-1)))*np.sin(np.radians(x[i]))
#         yyy=(8/(5*np.pi))*(
#             2*np.arccos((np.sqrt(2)+1)*xx/np.sqrt(2)) + np.arccos((np.sqrt(2)+1)*yy/np.sqrt(2))
#             )*(1-np.cos(np.radians(x[i])))
#         yi=zzz+xxx+yyy
#     y.append(yi)
    
# plt.plot(x,np.array(y), color = "red",label = "Mackenzi",linewidth=3)  
# plt.xlabel('Angle disorientation (deg)',fontsize=16)   
# plt.yticks(fontsize=15)
# plt.xticks(fontsize=15)
# plt.legend(fontsize=15) 
# plt.savefig("C:/Users/Zero/Documents/figures_for_paper/Disorientation_angle_Initial.svg",format="svg")
# # plt.show()grain_boundary_property_for_all_timesteps[0]       
# from svglib.svglib import svg2rlg
# from reportlab.graphics import renderPDF
# drawing = svg2rlg("C:/Users/Zero/Documents/figures_for_paper/Texture_Initial_after_equilibriation.svg")
# renderPDF.drawToFile(drawing, "C:/Users/Zero/Documents/PhD thesis/Curvature_effect/InitialMicrostructure/Texture_Initial_after_equilibriation.pdf")

'''
Normal=[]
for i in range(0,len(gId)):
    gbpi=grain_boundary_property_extended_for_all_timesteps[ts[i]]
    gbmesh=[j[1] for j in gbpi if j[0]==gId[i]][0]
    n=np.mean(gbmesh.vertex_normals,axis=0).tolist()
    Normal.append(n)
    

import math
from math import atan2, pi
def Quaternion2Euler(q):
        """
        Compute Euler angles from a Quaternion
        :param q: Quaternion
        :return: Euler angles (in degrees, Bunge convention)
        """
        P = 1
        (q0, q1, q2, q3) = q[0],q[1],q[2],q[3]
        q03 = q0 ** 2 + q3 ** 2
        q12 = q1 ** 2 + q2 ** 2
        chi = np.sqrt(q03 * q12)
        if chi == 0.:
            if q12 == 0.:
                phi_1 = atan2(-2 * P * q0 * q3, q0 ** 2 - q3 ** 2)
                Phi = 0.
            else:
                phi_1 = atan2(-2 * q1 * q2, q1 ** 2 - q2 ** 2)
                Phi = pi
            phi_2 = 0.
        else:
            phi_1 = atan2((q1 * q3 - P * q0 * q2) / chi,
                          (-P * q0 * q1 - q2 * q3) / chi)
            Phi = atan2(2 * chi, q03 - q12)
            phi_2 = atan2((P * q0 * q2 + q1 * q3) / chi,
                          (q2 * q3 - P * q0 * q1) / chi)+2*pi
        return [phi_1, Phi, phi_2]