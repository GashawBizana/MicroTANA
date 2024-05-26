import numpy as np
import trimesh
import copy


def Compute_change_in_time_Grain(M0,M1,G_id0,grain_mapping_data,cubMis,Method='ray'):
    
    
    def area_volume(m):
        mesh=copy.deepcopy(m)
        mesh.process()
        trimesh.repair.fill_holes(mesh)
        trimesh.repair.fix_normals(mesh, multibody=False)
        trimesh.repair.fix_inversion(mesh, multibody=False)
        
        return(mesh.area,mesh.volume)
    
    def initial_grain_id(grain_mapping_data, grain_id, time_step_n):
        return int(np.where(grain_mapping_data[time_step_n] == grain_id)[0][0])
    
    
    def grain_id_mapping_srl(grain_mapping_data, grain_id, time_step_n, step):
        index_n0 = np.where(grain_mapping_data[time_step_n] == grain_id)
        grain_id_m = np.nan if len(index_n0[0]) == 0 else grain_mapping_data[time_step_n + step][index_n0][0]
        return (np.nan,) if np.isnan(grain_id_m) else (int(grain_id_m),)
    
    def grains_recenter(grain_mesh_n_0, grain_mesh_m_1, edge_length):
        translation = np.rint((np.mean(grain_mesh_n_0.vertices, axis=0) - np.mean(grain_mesh_m_1.vertices, axis=0)) / edge_length) * edge_length
        grain_mesh_m_1.vertices += translation
        return grain_mesh_m_1

    def compute_distance(m1, m2, edge_length,method='ray'):
        m1=copy.deepcopy(m1)
        m2=copy.deepcopy(m2)
        m2 = grains_recenter(m1,m2, edge_length)
        m2.fix_normals()
        if method=='ray':
            
            ray_origins = m1.vertices
            
            inside = m2.ray.contains_points(ray_origins)
            sign = (inside.astype(int) * 2) - 1.0
            
            ray_directions = sign*m1.vertex_normals
            
            ray_idx=[i for i in range(len(ray_origins))]
            locations, index_ray, index_tri = m2.ray.intersects_location(
                ray_origins=ray_origins, ray_directions=ray_directions,multiple_hits=False
            )
            locations_all=np.array([ [np.nan,np.nan,np.nan]]*len(ray_idx))
            locations_all[index_ray]=locations
            signed_distances = np.linalg.norm(locations_all- ray_origins, axis=1)
            signed_distances *= sign
            
        
        elif method=='proximity':
            signed_distances = trimesh.proximity.signed_distance(m2, m1.vertices)
            
        computed_properties={'avg_Signed_distance':np.nanmean(signed_distances),'Signed_distances':signed_distances}
        
        return computed_properties
    
    
    G_id_initial0 = initial_grain_id(grain_mapping_data, G_id0[0], M0.timestep)
    detla_step=int(M1.timestep-M0.timestep)
    G_id1=grain_id_mapping_srl(grain_mapping_data, G_id0[0], M0.timestep, detla_step)
    
    
        
    
    
    if np.isnan(G_id1[0]):
        M0.grains[f'{G_id0}'].properties['Gid_1_exist']=False
        
        output={}
    
    else:
        
        M0.grains[f'{G_id0}'].properties['Gid_1_exist']=True
        
        a0,V0 =area_volume(M0.grains[f'{G_id0}'].properties['mesh'])
        a1,V1=area_volume(M1.grains[f'{G_id1}'].properties['mesh'])
        dV=V1-V0
        
        M0.grains[f'{G_id0}'].properties['dV']=dV
        
        area_mean_delta = 0.5 * (np.abs(a0) + np.abs(a1))
        
        M0.grains[f'{G_id0}'].properties['area_mean_delta']=area_mean_delta
        
        computed_properties=compute_distance(M0.grains[f'{G_id0}'].properties['mesh'], M1.grains[f'{G_id1}'].properties['mesh'], M1.EdgeLength,method=Method)
        
        avg_Signed_distance, Signed_distances=computed_properties['computed_properties'],computed_properties['Signed_distances']
        
        M0.grains[f'{G_id0}'].properties['avg_Signed_distance']=avg_Signed_distance
        M0.grains[f'{G_id0}'].properties['Signed_distances']=Signed_distances
        
        output={'Gid':G_id0,'Gid_initial':(G_id_initial0,),'dV':dV,'area_mean_delta':area_mean_delta,'avg_Signed_distance':avg_Signed_distance,'Signed_distances':Signed_distances}
    
    return output
    

    
def Compute_change_in_time_GrainBoundary(M0,M1,GB_id0,grain_mapping_data,Use_grain=True):
    
    
    def area_compute(m):
        mesh=copy.deepcopy(m)
        mesh.process()
        trimesh.repair.fix_normals(mesh, multibody=False)
        trimesh.repair.fix_inversion(mesh, multibody=False)
        
        return(mesh.area)
    
    
    def grainboundary_recenter(grainboundary_mesh_n_0, grainboundary_mesh_m_1, edge_length):
        translation = np.rint((np.mean(grainboundary_mesh_n_0.vertices, axis=0) - np.mean(grainboundary_mesh_m_1.vertices, axis=0)) / edge_length) * edge_length
        grainboundary_mesh_m_1.vertices += translation
        return grainboundary_mesh_m_1
    
    def distance_from_Grain(M0,GBid):
        
        triangle01,triangle10=M0.grain_boundaries[f'{GBid}'].properties['triangle_id']
        
        Signed_distances0=M0.grains[f'{(GBid[0],)}'].properties['Signed_distances']
        Signed_distances1=M0.grains[f'{(GBid[1],)}'].properties['Signed_distances']
        
        Signed_distances01=1/3*(Signed_distances0[triangle01[0,:]]+Signed_distances0[triangle01[1,:]] +Signed_distances0[triangle01[2,:]])
        avg_Signed_distance01=np.nanmean(Signed_distances01)
        
        Signed_distances10=1/3*(Signed_distances1[triangle10[0,:]]+Signed_distances1[triangle10[1,:]] +Signed_distances1[triangle10[2,:]])
        avg_Signed_distance10=np.nanmean(Signed_distances10)
        
        computed_properties={'Signed_distances':(Signed_distances01,Signed_distances10),'avg_Signed_distance':(avg_Signed_distance01,avg_Signed_distance10)}
        
        return computed_properties
        
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
            
            grain_boundary_forming_grains_m=(m0,m1)
        else:
            
            grain_boundary_forming_grains_m=(int(m0),int(m1))
        
        return(grain_boundary_forming_grains_m)    
    
    def Initial_Grain_Id(grain_mapping_data, grain_boundary_forming_grains_n,time_step_n):
        n0=grain_boundary_forming_grains_n[0]
        n1=grain_boundary_forming_grains_n[1]
        grain_mapping_data_time_step_n=grain_mapping_data[time_step_n]
        index_n0=np.where(grain_mapping_data_time_step_n==n0)
        index_n1=np.where(grain_mapping_data_time_step_n==n1)
        
        return (int(index_n0[0]),int(index_n1[0]))
    
    def compute_distance_GB(m1,m2, edge_length):
        
        m1=copy.deepcopy(m1)
        m2=copy.deepcopy(m2)
        
        m2 = grainboundary_recenter(m1, m2, edge_length)
        
        ray_origins = m1.vertices
        
        inside = m2.ray.contains_points(ray_origins)
        sign = (inside.astype(int) * 2) - 1.0
        
        ray_directions = sign*m1.vertex_normals
        
        ray_idx=[i for i in range(len(ray_origins))]
        locations, index_ray, index_tri = m2.ray.intersects_location(
            ray_origins=ray_origins, ray_directions=ray_directions,multiple_hits=False
        )
        locations_all=np.array([ [np.nan,np.nan,np.nan]]*len(ray_idx))
        locations_all[index_ray]=locations
        signed_distances = np.linalg.norm(locations_all- ray_origins, axis=1)
        signed_distances *= sign
        
        return (signed_distances)
    
    
    #######################################
    detla_step=int(M1.timestep-M0.timestep)
    
    GB_id1=GrainIdMapper(grain_mapping_data, GB_id0,M0.timestep,detla_step)
    GB_id_initial0=Initial_Grain_Id(grain_mapping_data, GB_id0,M0.timestep)
    
        
        

        
    if np.isnan(GB_id1[0]) or np.isnan(GB_id1[0]):
        M0.grainboundaries[f'{GB_id0}'].properties['GBid_1_exist']=False
        
        output={}
    
    else:
        
        M0.grainboundarie[f'{GB_id0}'].properties['GBid_1_exist']=True
            
        if Use_grain==True:
            computed_properties=distance_from_Grain(M0,GB_id0)
            
        if Use_grain==False:
            
            
            mesh_0_0,mesh_0_1=M0.grains[f'{GB_id0}'].properties['mesh']
            
            
            
            if f'{GB_id1}' not in M1.grains.keys:
                GB_id1=(GB_id1[1],GB_id1[0])
                
                a_0_0,a_0_1=area_compute(M0.grainboundarie[f'{GB_id0}'].properties['mesh'][0]),area_compute(M0.grainboundarie[f'{GB_id0}'].properties['mesh'][1])
                a_1_0,a_1_1=area_compute(M1.grainboundarie[f'{GB_id1}'].properties['mesh'][0]),area_compute(M1.grainboundarie[f'{GB_id1}'].properties['mesh'][1])
            
            area_mean_delta_00 = 0.5 * (np.abs(a_0_0) + np.abs(a_1_0))
            
            area_mean_delta_11= 0.5 * (np.abs(a_0_1) + np.abs(a_1_1))
            
            
            M0.grainboundarie[f'{GB_id0}'].properties['area_mean_delta']=(area_mean_delta_00,area_mean_delta_11)
            
            mesh_1_0,mesh_1_1=M1.grains[f'{GB_id1}'].properties['mesh']
            
            Signed_distances00=compute_distance_GB(mesh_0_0,mesh_1_0, M1.EdgeLength)
            
            Signed_distances11=compute_distance_GB(mesh_0_1,mesh_1_1, M1.EdgeLength)
            
            computed_properties={'avg_Signed_distance':(np.nanmean(Signed_distances00),np.nanmean(Signed_distances11)),'Signed_distances':(Signed_distances00,Signed_distances11)}
            
            M0.grainboundarie[f'{GB_id0}'].properties['avg_Signed_distance']=computed_properties['avg_Signed_distance']
            M0.grainboundarie[f'{GB_id0}'].properties['Signed_distances']=computed_properties['Signed_distances']
            
            output={'GBid':GB_id0,'GBid_initial':GB_id_initial0,'area_mean_delta':(area_mean_delta_00,area_mean_delta_11),'avg_Signed_distance':computed_properties['avg_Signed_distance'],'Signed_distances':computed_properties['Signed_distances']}
        
        return output