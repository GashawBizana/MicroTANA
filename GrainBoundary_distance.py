import pyvista as pv
import numpy as np
import time

import matplotlib as mpl
import matplotlib.pylab as plt


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

def distanace_between_current_and_next(grain_boundary_step_i, grain_boundary_step_i_plus_1,EdgeLength_step_i_plus_1):
    
    grain_boundary_step_i_cm=np.mean(grain_boundary_step_i.vertices,axis=0)
    grain_boundary_step_i_plus_1_cm=np.mean(grain_boundary_step_i_plus_1.vertices,axis=0)
    
    diff_center=np.rint((grain_boundary_step_i_cm-grain_boundary_step_i_plus_1_cm)/EdgeLength_step_i_plus_1)    
    periodic_translation_factor=diff_center*EdgeLength_step_i_plus_1
    grain_boundary_step_i_plus_1.vertices= grain_boundary_step_i_plus_1.vertices + periodic_translation_factor
    
    # cloud_grain_boundary_step_i= pv.PolyData(grain_boundary_step_i)
    # cloud_grain_boundary_step_i_plus_1=pv.PolyData(grain_boundary_step_i_plus_1)
    
    
    
    # h0= cloud_grain_boundary_step_i.delaunay_2d()
    # h1 = cloud_grain_boundary_step_i_plus_1.delaunay_2d()
    
    h0=pv.wrap(grain_boundary_step_i)
    h1=pv.wrap(grain_boundary_step_i_plus_1)
    '''
    distance=np.sqrt(np.sum((np.array(h0.center) - np.array(h1.center)) ** 2))
    ''' 
    '''
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
    # mean=np.nanmean(h0n["distances"])
    # std=np.nanstd(h0n["distances"])
    # distance= np.nanmean(h0n["distances"][h0n["distances"]  <= 100])
    distance= np.nanmean(h0n["distances"])
    # distance=np.nanmean( np.where(h0n["distances"] <= 2, np.nan, h0n["distances"]))
    
    
    '''
    
    closest_cells, closest_points = h1.find_closest_cell(h0.points, return_closest_point=True)
    d_exact = np.linalg.norm(h0.points - closest_points, axis=1)
    h0["distances"] = d_exact
    distance=np.mean(d_exact)
    
    
    
    curv0 = h0.curvature()
    curv1=h1.curvature()
    average_mean_curvature =np.mean(np.abs(curv0))
    
    dA=np.abs(h1.area-h0.area)
    
    
    return(distance,average_mean_curvature,dA,[h0,h1])



def GrainBoundary_distance_computation(EdgeLength_for_all_timesteps,grain_boundary_property_for_all_timesteps,GrainBoundary_length_of_triple_line_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps,GrainBoundary_number_of_triple_line_for_all_timesteps,ii,grain_mapping_data,step):

    Grain_Boundary_Velocity_Data_for_all_timesteps=[]
    
    for k in range(0,len(grain_boundary_property_for_all_timesteps)-step):
        
        EdgeLength_step_k= EdgeLength_for_all_timesteps[k]
        EdgeLength_step_k_plus_1=EdgeLength_for_all_timesteps[k+step]
        grain_boundary_property_step_k=grain_boundary_property_for_all_timesteps[k]
        grain_boundary_property_step_k_plus_1=grain_boundary_property_for_all_timesteps[k+step]
        grain_boundary_forming_step_k=[j[0] for j in grain_boundary_property_step_k]
        
    
        
        grain_boundary_mean_curvature_step_k=[j[3] for j in grain_boundary_property_step_k]
        grain_boundary_mean_curvature_step_k_plus_1=[j[3] for j in grain_boundary_property_step_k_plus_1]
    
        grain_boundary_total_mean_curvature_step_k=[j[1] for j in grain_boundary_property_step_k]
        grain_boundary_total_mean_curvature_step_k_plus_1=[j[1] for j in grain_boundary_property_step_k_plus_1]
        
        grain_boundary_total_surface_area_step_k=[j[2] for j in grain_boundary_property_step_k]
        grain_boundary_total_surface_area_step_k_plus_1=[j[2] for j in grain_boundary_property_step_k_plus_1]
        
        grain_boundary_mesh_step_k=[j[4] for j in grain_boundary_property_step_k]
        grain_boundary_mesh_step_k_plus_1=[j[4] for j in grain_boundary_property_step_k_plus_1]
        
        GrainBoundary_length_of_triple_line_step_k=GrainBoundary_length_of_triple_line_for_all_timesteps[k]
        GrainBoundary_length_of_triple_line_step_k_plus_1=GrainBoundary_length_of_triple_line_for_all_timesteps[k+step]
        GrainBoundary_misorientation_step_k=GrainBoundary_misorientation_for_all_timesteps[k]
        GrainBoundary_misorientation_step_k_plus_1=GrainBoundary_misorientation_for_all_timesteps[k+step]
        GrainBoundary_number_of_triple_line_step_k=GrainBoundary_number_of_triple_line_for_all_timesteps[k]
        GrainBoundary_number_of_triple_line_step_k_plus_1=GrainBoundary_number_of_triple_line_for_all_timesteps[k+step]
        
        Grain_Boundary_Velocity_Data_k=[]
        for i in range(0,len(grain_boundary_property_step_k)):
            
            grain_id_k_i=grain_boundary_forming_step_k[i]
            grain_id_k_i_initial=Initial_Grain_Id(grain_mapping_data, grain_id_k_i ,k)
            
            grain_id_k_plus_1_i=GrainIdMapper(grain_mapping_data, grain_id_k_i ,k,step)
            
            grain_boudary_position_step_k_i=grain_boundary_mesh_step_k[i]
            
            grain_boundary_mean_curvature_step_k_i=grain_boundary_mean_curvature_step_k[i]
            
            grain_boundary_total_mean_curvature_step_k_i=grain_boundary_total_mean_curvature_step_k[i]
            
            grain_boundary_total_surface_area_step_k_i=grain_boundary_total_surface_area_step_k[i]
            
            GrainBoundary_length_of_triple_line_step_k_i=[j[1] for j in GrainBoundary_length_of_triple_line_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            GrainBoundary_misorientation_step_k_i= [j[1] for j in GrainBoundary_misorientation_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            GrainBoundary_number_of_triple_line_step_k_i= [j[1] for j in GrainBoundary_number_of_triple_line_step_k if (grain_id_k_i[0] in j[0] and grain_id_k_i[1] in j[0]) ][0]
            
            
            if (np.isnan(grain_id_k_plus_1_i).any()) ==False:
                grain_boudary_position_k_plus_1_i=[j[4] for j in grain_boundary_property_step_k_plus_1 if j[0]==grain_id_k_plus_1_i]
                
                if (not grain_boudary_position_k_plus_1_i )==False:
                    
                    
                    
                    grain_boudary_position_k_plus_1_i=grain_boudary_position_k_plus_1_i[0]
                    
                    distance_step_k_j,average_mean_curvature_pyvista_step_k_j,dA_step_k_j,h0h1=distanace_between_current_and_next(grain_boudary_position_step_k_i, grain_boudary_position_k_plus_1_i,EdgeLength_step_k_plus_1)
                    
                    Grain_Boundary_Velocity_Data_k.append([grain_id_k_i, distance_step_k_j,average_mean_curvature_pyvista_step_k_j,grain_boundary_mean_curvature_step_k_i,grain_boundary_total_mean_curvature_step_k_i, GrainBoundary_misorientation_step_k_i,GrainBoundary_number_of_triple_line_step_k_i,GrainBoundary_length_of_triple_line_step_k_i,grain_boundary_total_surface_area_step_k_i,dA_step_k_j,h0h1,grain_id_k_plus_1_i ])
        
        
        Grain_Boundary_Velocity_Data_for_all_timesteps.append(Grain_Boundary_Velocity_Data_k)
        
    return(Grain_Boundary_Velocity_Data_for_all_timesteps)
step=1
t = time.time()
Grain_Boundary_Velocity_Data_for_all_timesteps=GrainBoundary_distance_computation(EdgeLength_for_all_timesteps,grain_boundary_property_for_all_timesteps,GrainBoundary_length_of_triple_line_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps,GrainBoundary_number_of_triple_line_for_all_timesteps,ii,grain_mapping_data,step)

'''
Output of Grain_Boundary_Velocity_Data_for_all_timesteps is as follow:
    
    0 -grain_id, 1- distance ,2 - average_mean_curvature_pyvista , 3 -grain_boundary_mean_curvature_trimesh, 4- grain_boundary_total_mean_curvature_trimesh, 5 - GrainBoundary_misorientation ,6 - GrainBoundary_number_of_triple_line,7 -  GrainBoundary_length_of_triple_line , 8 - grain_boundary_total_surface_area, 9- dA
    

'''

print (" Velocity Analysia took: ", (time.time()-t))

curvlimit=[0,1]
dislimit=[0,100]
mislimit=[0,66]
arealimit=[]
grain1=[0,100]
grain2=[0,0]

LS=0
LE=1#len(Grain_Boundary_Velocity_Data_for_all_timesteps)

dis,curv_pyvista,curv,mis,tot_curv,tot_area,TL_n,TL_l, dA,h0h1=[],[],[],[],[],[],[],[],[],[]

for i in range(LS,LE):
    
    
    dis0=[k[1] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1]) ]
    
    
    curv_pyvista0=[k[2] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]
    
    
    curv0=[k[3] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]
    
    tot_curv0=[k[4] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]  
    
    mis0=[k[5] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]
    
    TL_n0=[k[6] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]
    
    TL_l0=[k[7] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]
    
    tot_area0=[k[8] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]
    
    dA0=[k[9] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]
    h0h10=[k[10] for k in Grain_Boundary_Velocity_Data_for_all_timesteps[i] if (np.isnan( k[1]).any()) ==False and dislimit[0] < k[1] < dislimit[1] and curvlimit[0] < k[3] <curvlimit[1] 
and mislimit[0] < k[5] < mislimit[1] and (grain1[0]<=k[0][0]<=grain1[1] or grain1[0]<=k[0][0]<=grain1[1])]
    
    dis=dis+dis0
    curv_pyvista=curv_pyvista+curv_pyvista0
    curv=curv+curv0
    tot_curv=tot_curv+tot_curv0
    mis=mis+mis0
    TL_n=TL_n+TL_n0
    TL_l=TL_l+TL_l0
    tot_area=tot_area+tot_area0
    dA=dA+dA0
    h0h1=h0h1+h0h10


view=np.array(mis)
xa=np.array(curv_pyvista)
ya=np.array(dis)*0.1

cmap=plt.cm.jet
cmaplist=[cmap(i) for i in range(cmap.N)]
cmaplist[0]=(0.5,0.5,0.5,1.0)
cmap=mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist,cmap.N)

# bounds=np.linspace(np.rint(min(mis)),np.rint(max(mis)),2*int(np.rint(max(mis))-np.rint(min(mis))))
bounds=np.linspace(min(view),max(view),256)
norm=mpl.colors.BoundaryNorm(bounds, cmap.N)
coef = np.polyfit(xa,ya,1)
poly1d_fn = np.poly1d(coef)
xnew=np.linspace(min(xa),max(xa),10)
fig = plt.figure()
ax = fig.add_subplot()

sc=ax.scatter(xa,ya,c=view,cmap=cmap,norm=norm)

ax.set_xlabel('curvature (1/$\AA$)')
# ax.set_ylabel('phi')
ax.set_ylabel('velocity ($\AA$/time_step)')
cb=plt.colorbar(sc)
cb.set_label('misorientation (deg)')
plt.plot(xnew,poly1d_fn(xnew))
plt.show()