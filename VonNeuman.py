import numpy as np

import time

import matplotlib.pylab as plt
import matplotlib as mpl
import trimesh


def Initial_Grain_Id(grain_mapping_data, grain_id,time_step_n):
    grain_mapping_data_time_step_n=grain_mapping_data[time_step_n]
    index_n0=np.where(grain_mapping_data_time_step_n==grain_id)
    return(int(index_n0[0]))

def GrainId_mapping_srl(grain_mapping_data, grain_id,time_step_n,step):

    
    grain_mapping_data_time_step_n=grain_mapping_data[time_step_n]
    index_n0=np.where(grain_mapping_data_time_step_n==grain_id)
    
    grain_mapping_data_time_step_m=grain_mapping_data[time_step_n+step]
    
    if len(index_n0[0])==0:
        m0=np.nan
    else:
        m0=(grain_mapping_data_time_step_m[index_n0])[0]
    
    if np.isnan(m0):
        
        grain_id_m=m0
    else:
        
        grain_id_m=int(m0)
        
    return(grain_id_m)

def GrainsRecenter(grain_mesh_n_0,grain_mesh_m_1,EdgeLength_step_i_plus_1):
    
    grain_mesh_n_0_cm=np.mean(grain_mesh_n_0.vertices,axis=0)
    grain_mesh_m_1_cm=np.mean(grain_mesh_m_1.vertices,axis=0)
    
    diff_center=np.rint((grain_mesh_n_0_cm-grain_mesh_m_1_cm)/EdgeLength_step_i_plus_1)    
    periodic_translation_factor=diff_center*EdgeLength_step_i_plus_1
    grain_mesh_m_1.vertices= grain_mesh_m_1.vertices + periodic_translation_factor
    
    return(grain_mesh_m_1)


def distance(m1,m2,EdgeLength_step_k_plus_1):
    
    m2=GrainsRecenter(m1,m2,EdgeLength_step_k_plus_1)
    # t=trimesh.proximity.longest_ray(m2, m1.vertices, -1*m1.vertex_normals)
    # t1=trimesh.proximity.longest_ray(m2, m1.vertices, m1.vertex_normals)
    tmc=(1/(2*np.pi))*(m1.integral_mean_curvature)
    # tt=[]
    # for i in range(0,len(t)):
    #     dt=t[i]
    #     dt1=t1[i]
    #     if dt==np.inf and dt1==np.inf:
    #         tt.append(np.nan)
    #     elif dt <=dt1:
    #         tt.append(-dt)
    #     else:
    #         tt.append(dt1)
    tt=trimesh.proximity.signed_distance(m2, m1.vertices)
    area_mean=0.5*(np.abs(m1.area)+np.abs(m2.area))
    # tt=np.ones(len(m1.vertices))
    dis=np.nanmean(tt)
    
    return(dis,tmc, tt,area_mean)

def Grain_VonNeuman_Data (grain_property_for_all_timesteps,grain_total_curvature_surface_area_volume_for_all_timesteps,grain_total_curvature_surface_area_volume_from_face_for_all_timesteps,Grain_average_misorientation_for_all_timesteps,Grain_number_of_grain_boundaries_for_all_timesteps,Grain_number_of_triple_lines_for_all_timesteps,Grain_length_of_triple_lines_for_all_timesteps,grain_mapping_data,step,EdgeLength_for_all_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps):
    
    Grain_VonNeuman_Data_for_timesteps=[]
    # len(grain_property_for_all_timesteps)-step
    for k in range(0,len(grain_property_for_all_timesteps)-step):
        print(k)
        
        GrainBoundary_total_curvature_mean_curvature_surface_area_step_k=GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps[k]
        GrainBoundary_misorientation_step_k=GrainBoundary_misorientation_for_all_timesteps[k]
        
        EdgeLength_step_k= EdgeLength_for_all_timesteps[k]
        EdgeLength_step_k_plus_1=EdgeLength_for_all_timesteps[k+step]
        
        grain_property_step_k=grain_property_for_all_timesteps[k]
        
        grain_total_curvature_step_k=[j[1] for j in grain_total_curvature_surface_area_volume_for_all_timesteps[k]]
        grain_total_surface_area_step_k=[j[2] for j in grain_total_curvature_surface_area_volume_for_all_timesteps[k]]
        grain_total_volume_step_k =[j[3] for j in grain_total_curvature_surface_area_volume_for_all_timesteps[k]]
        
        grain_total_curvature_face_step_k=[j[1] for j in grain_total_curvature_surface_area_volume_from_face_for_all_timesteps[k]]
        
        grain_total_surface_area_face_step_k=[j[2] for j in grain_total_curvature_surface_area_volume_from_face_for_all_timesteps[k]]
        
        
        
        grain_total_volume_face_step_k =[j[3] for j in grain_total_curvature_surface_area_volume_from_face_for_all_timesteps[k]]
        grain_total_curvature_face_positive_step_k =[j[4] for j in grain_total_curvature_surface_area_volume_from_face_for_all_timesteps[k]]
        
        Grain_average_misorientation_step_k=Grain_average_misorientation_for_all_timesteps[k]
        Grain_number_of_grain_boundaries_step_k=Grain_number_of_grain_boundaries_for_all_timesteps[k]
        Grain_number_of_triple_lines_step_k=Grain_number_of_triple_lines_for_all_timesteps[k]
        Grain_length_of_triple_lines_step_k=Grain_length_of_triple_lines_for_all_timesteps[k]
        
        
        Grain_VonNeuman_Data_step_k=[]
        for j in range(0,len(grain_property_step_k)):
            
            
            
            grain_id_k_j=grain_property_step_k[j][0]
            grain_id_k_j_initial=Initial_Grain_Id(grain_mapping_data,grain_id_k_j,k)
            grain_mesh_k_j=grain_property_step_k[j][1]
            
            grain_total_curvature_k_j=grain_total_curvature_step_k[j]
            grain_total_surface_area_k_j=grain_total_surface_area_step_k[j]
            grain_total_volume_k_j= grain_total_volume_step_k[j]
            
            grain_total_curvature_face_k_j=grain_total_curvature_face_step_k[j]
            grain_total_surface_area_face_k_j=grain_total_surface_area_face_step_k[j]
            grain_total_volume_face_k_j= grain_total_volume_face_step_k[j]
            grain_total_curvature_face_positive_k_j=grain_total_curvature_face_positive_step_k[j]
            
            Grain_average_misorientation_k_j=Grain_average_misorientation_step_k[j][1]
            Grain_number_of_grain_boundaries_k_j= Grain_number_of_grain_boundaries_step_k[j][1]
            
            Extended_neighbour=[]
            for IDX in range(0,len(Grain_number_of_grain_boundaries_step_k[j][2])):
                gbId=Grain_number_of_grain_boundaries_step_k[j][2][IDX]
                gn=[gni for gni in gbId if gni !=grain_id_k_j][0]
                en=[gbni[1] for gbni in Grain_number_of_grain_boundaries_step_k if gbni[0][0] == gn]
                if not en:
                    Extended_neighbour.append(np.nan)
                else:
                    Extended_neighbour.append(en[0])
            Average_gbn_neighbour=np.nanmean(Extended_neighbour)
            
            area_of_gb_associated_with_grain=[j[3] for j in GrainBoundary_total_curvature_mean_curvature_surface_area_step_k if (grain_id_k_j == j[0][0])]
            
            misorientation_of_gb_associated_with_grain=[j[1] for j in GrainBoundary_misorientation_step_k if (grain_id_k_j in j[0])]
            wighted_area_of_gb_associated_with_grain=np.array(area_of_gb_associated_with_grain)/np.nansum(np.array(area_of_gb_associated_with_grain))
            
            area_wighted_average_misorientation_of_gb_associated_with_grain=np.array(misorientation_of_gb_associated_with_grain)*wighted_area_of_gb_associated_with_grain
            mis_weighted=np.nansum(area_wighted_average_misorientation_of_gb_associated_with_grain)
            
            Grain_number_of_triple_lines_k_j= Grain_number_of_triple_lines_step_k[j][1]
            Grain_length_of_triple_lines_k_j= Grain_length_of_triple_lines_step_k[j][1]
            
            
            grain_id_k_plus_1_j=GrainId_mapping_srl(grain_mapping_data,grain_id_k_j,k,step)
            grain_mesh_step_k_plus_1= [i for i in grain_property_for_all_timesteps[k+step]]
            grain_mesh_k_plus_1_j=[i[1] for i in grain_mesh_step_k_plus_1 if i[0]==grain_id_k_plus_1_j]
            
            grain_total_surface_area_face_k_plus_1_j=[i[2] for i in grain_total_curvature_surface_area_volume_from_face_for_all_timesteps[k+step] if i[0]==grain_id_k_plus_1_j]
            
            grain_total_volume_k_plus_1_j= [i[3] for i in grain_total_curvature_surface_area_volume_for_all_timesteps[k+step] if i[0]==grain_id_k_plus_1_j]
        
            
            
            if (not grain_mesh_k_plus_1_j )==False:
                
                dv=(-np.abs(grain_total_volume_k_j)+np.abs(grain_total_volume_k_plus_1_j))[0]
                dis,tmc,tt,area_mean=distance(grain_mesh_k_j,grain_mesh_k_plus_1_j[0],EdgeLength_step_k_plus_1)
                area_mean_face=0.5*(np.abs(grain_total_surface_area_face_k_plus_1_j[0])+np.abs(grain_total_surface_area_face_k_j))
            # elif (not grain_mesh_k_plus_1_j )==True:
                
            #     dv=-np.abs(grain_total_volume_k_j)
            #     dis,tmc,tt=(-np.abs(grain_total_volume_k_j)+np.abs(0))/grain_mesh_k_j.area,(1/(2*np.pi))*(grain_mesh_k_j.integral_mean_curvature),np.ones(len(grain_mesh_k_j.vertices))*(-np.abs(grain_total_volume_k_j)+np.abs(0))/grain_mesh_k_j.area
                
                Grain_VonNeuman_Data_step_k.append([grain_id_k_j_initial,grain_total_curvature_k_j, dv, Grain_length_of_triple_lines_k_j, Grain_number_of_triple_lines_k_j,Grain_number_of_grain_boundaries_k_j,Grain_average_misorientation_k_j,grain_total_volume_k_j, grain_total_surface_area_k_j,grain_total_curvature_face_k_j,grain_total_surface_area_face_k_j, grain_total_volume_face_k_j,dis,tmc,tt,grain_total_curvature_face_positive_k_j,area_mean,area_mean_face,Average_gbn_neighbour,mis_weighted])
            
        
    
        Grain_VonNeuman_Data_for_timesteps.append(Grain_VonNeuman_Data_step_k)
    
    
    return(Grain_VonNeuman_Data_for_timesteps)
                
    
    '''
    Output of Grain_VonNeuman_Data_for_timesteps is as follow:
        
        0 -grain_id, 1- grain_total_curvature ,2 - change in volume(dV) , 3 -Grain_length_of_triple_lines, 4- Grain_number_of_triple_lines, 5 - Grain_number_of_grain_boundaries ,6 - Grain_average_misorientation,7 -  grain_total_volume , 8 - grain_total_surface_area
        

    '''
step=1
t = time.time()
Grain_VonNeuman_Data_for_timesteps=Grain_VonNeuman_Data (grain_property_for_all_timesteps,grain_total_curvature_surface_area_volume_for_all_timesteps,grain_total_curvature_surface_area_volume_from_face_for_all_timesteps,Grain_average_misorientation_for_all_timesteps,Grain_number_of_grain_boundaries_for_all_timesteps,Grain_number_of_triple_lines_for_all_timesteps,Grain_length_of_triple_lines_for_all_timesteps,grain_mapping_data,step,EdgeLength_for_all_timesteps,GrainBoundary_total_curvature_mean_curvature_surface_area_for_all_timesteps,GrainBoundary_misorientation_for_all_timesteps)

print (" Grain volume change analysis took : ", (time.time()-t))
# Grain_VonNeuman_Data_for_timesteps=Grain_VonNeuman_Data_for_timesteps_with_velocity

G_tc, G_dv, G_TL_l, G_TL_n, G_GB_n, G_mis, G_v, G_a, G_id, G_tcf,G_step,G_dis,G_tmc,G_af,G_tcf_p,G_a_mean,G_a_meanf,G_GB_n_ext,G_mis_w= [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]

g=[0,4]
tc=[-10000000000,10000000000]
dis=[-100000,10000]
tln=[0,2000]
gbn=[0,100]
LS=0

LE=len(Grain_VonNeuman_Data_for_timesteps)
for i in range(0,LE):
    
    G_id0=[j[0] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]] 
    G_tc0=[j[1] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_dv0=[j[2] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_TL_l0=[j[3] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_TL_n0=[j[4] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_GB_n0=[j[5] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_mis0=[j[6] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_v0=[j[7] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_a0=[j[8] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_tcf0=[j[9] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_step0=[i for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_dis0=[j[12] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_tmc0=[j[13] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_af0=[j[10] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_tcf_p0=[j[15] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_a_mean0=[j[16] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_a_meanf0=[j[17] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_GB_n_ext0=[j[18] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    G_mis_w0=[j[19] for j in Grain_VonNeuman_Data_for_timesteps[i] if g[0]<=j[0] <= g[1] and dis[0]<=j[12] <= dis[1] and  tc[0] <= j[1] <= tc[1] and tln[0] <= j[4] <= tln[1]  and gbn[0] <= j[5] <= gbn[1]]
    
    G_id=G_id+G_id0
    G_tc=G_tc+G_tc0
    G_dv=G_dv+G_dv0
    G_TL_l=G_TL_l+G_TL_l0
    G_TL_n=G_TL_n+G_TL_n0
    G_GB_n=G_GB_n+G_GB_n0
    G_mis=G_mis+G_mis0
    G_v=G_v+G_v0
    G_a=G_a+G_a0
    G_tcf=G_tcf+G_tcf0
    G_step=G_step+G_step0
    G_dis=G_dis+G_dis0
    G_tmc=G_tmc+G_tmc0
    G_af=G_af+G_af0
    G_tcf_p=G_tcf_p+G_tcf_p0
    G_a_mean=G_a_mean+G_a_mean0
    G_a_meanf=G_a_meanf+G_a_meanf0
    G_GB_n_ext=G_GB_n_ext+G_GB_n_ext0
    G_mis_w=G_mis_w+G_mis_w0
view=np.array(G_mis_w)
dt=1/(step*0.01)
xa1=np.array(G_tc)
xa2=np.array(G_TL_l)
xa=-2*np.pi*(xa1-xa2/6)*0.1
xa3=np.array(G_v)**(1/3)*0.1
ya=-xa/(np.array(G_a) *0.01)

# critR=[]
# for k in range(0,len(xa)):
#         critR.append(np.mean(np.array([(9/8)*(j*0.001)**(1/3) for idx,j in enumerate(G_v) if G_step[idx]==G_step[k]])))
    

# ya=np.array(G_tcf)*0.1/(np.array(G_af) *0.01)
# view=xa
# ya=np.array(G_dv)*0.001/(np.array(G_a) *0.01)*np.power(np.array(G_v)*0.001, -1/3)
# ya=np.array(G_dis)*0.1*dt
# ya=np.power(np.array(G_v)*0.001, -1/3)*np.array(G_dv)*dt*0.001
# xa=10*np.array(G_tcf_p)/np.array(G_af)

# ya=np.array(G_dv)*dt*0.001*np.power(np.array(G_v)*0.001, -2/3)
# xa=np.power(np.array(G_v)*0.001, -1/3)
# ya=np.array(G_dis)*0.1*dt
# xa=10*np.array(G_tc-xa2/6)/np.array(G_a)

# xa=-(np.array(G_tc)-xa2/6)
# ya=np.array(G_dv)*0.001*np.power(np.array(G_v)*0.001, -2/3)*(1/(4*np.pi))
# xa=-10*np.array(G_tc-xa2/6)/np.array(G_a)
# xa=np.power(np.array(G_v)*0.001,-1/3)

# xa=-(np.array(G_tcf))*0.1/(np.array(G_af)*0.01)
# xa=np.array(G_TL_l)/6
# plt.rcParams['text.usetex'] = True
# xa=np.power(np.array(G_v)*0.001, 1/3)/(np.mean(np.power(np.array(G_v)*0.001, 1/3))*9/8)
# xa=np.array(G_GB_n)
# ya=np.power(np.array(G_v)*0.001,-1/3)*np.array(G_dv)*dt*0.001/(4*np.pi)
# xa=np.array(G_GB_n)
# ya=np.array(G_dv)*dt*0.001
# xa=np.array(G_tc)/np.array(G_a)
# ya=np.array(G_dv)*np.power(np.array(G_v)*0.001,-2/3)*0.001*dt*(3*0.25/np.pi)*dt
# ya=np.array(G_dis)*0.1*dt
xa=-2*np.pi*np.array(G_tcf)*0.1
# ya=np.array(G_dv)*0.001*dt
#ya=-2*np.pi*(xa1-xa2/6)*0.1
# ya=np.power(np.array(G_v)*0.001, 1/3)/np.array(critR)
#xa=np.power(np.array(G_v)*0.001, -1/3)
# xa=np.array(G_step)*0.01
# ya=(xa1)*0.1
#ya=(xa2/6)*0.1
# ya=np.array(G_dis)*0.1*dt
#xa=np.array(G_GB_n)-np.array(G_GB_n_ext)
ya=np.array(G_dv)*0.001*dt
# xa=np.array(G_step)
# xa=-2*np.pi*np.array(G_tcf)*0.1/(np.array(G_af)*0.01)
# ya=np.array(G_dv)*0.001*dt/(np.array(G_a_meanf)*0.01)/(0.01)
# ya=np.array(G_dis)*0.1*dt
# plt.rcParams['text.usetex'] = False
cmap=plt.cm.jet
cmaplist=[cmap(i) for i in range(cmap.N)]
cmaplist[0]=(0.5,0.5,0.5,1.0)
cmap=mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist,cmap.N)

# bounds=np.linspace(np.rint(min(mis)),np.rint(max(mis)),2*int(np.rint(max(mis))-np.rint(min(mis))))
bounds=np.linspace(min(view),max(view))
norm=mpl.colors.BoundaryNorm(bounds, cmap.N)
coef = np.polyfit(xa,ya,1)
poly1d_fn = np.poly1d(coef)
xnew=np.linspace(min(xa),max(xa),38)
fig = plt.figure()
ax = fig.add_subplot()


# sc=ax.scatter(xa,(1/6) *xa2*0.1,c=view,cmap=cmap,norm=norm,s=25,marker='+',label=r'$\frac{1}{6}L_{t}$')
# sc=ax.scatter(xa,xa1*0.1,c=view,cmap=cmap,norm=norm,s=25,marker='v',label=r'$\mathcal{L}$')
# ax.legend(loc="upper right")
# ax.scatter(xa3,xa2/6,c=view,cmap=cmap,norm=norm)
#ax.scatter(xa,ya1)
sc=ax.scatter(xa,ya,c=view,cmap=cmap,norm=norm,s=50,marker='+')

# ax.set_xlabel('M/6-L ($nm$)')
# ax.set_ylabel(r'$2\pi\left(\frac{1}{6}L_{t} - \mathcal{L}\right)(nm)$',fontweight='bold',size=12)
ax.set_ylabel(r'$\frac{1}{6}L_{t}(nm)$',fontweight='bold',size=12)
# ax.set_ylabel(r'$\mathcal{L}(nm)$',fontweight='bold',size=12)

# ax.set_xlabel(r'Average Grain Curvature($\frac{1}{nm}$)',fontweight='bold',size=12)
# ax.set_ylabel(r'Velocity($\frac{nm}{ns}$)',fontweight='bold',size=12)

ax.set_xlabel(r'$V^{-1/3} (nm)^{-1/3}$',fontweight='bold',size=12)
# ax.set_xlabel(r'$-2\pi\mathcal{L}_{Face}(nm)$',fontweight='bold',size=12)
# ax.set_ylabel(r'$\frac{R}{R_{cr}}$',fontweight='bold',size=12)
# ax.set_xlabel(r'$Time (ns)}$',fontweight='bold',size=12)
# ax.set_ylabel(r'$\frac{dV}{dt} (\frac{nm^3}{ns})$',fontweight='bold',size=12)
# ax.set_xlabel(r'$\frac{dV}{dt}A_{mean}^{-1} (nm/t_{step})$',fontweight='bold',size=12)
# ax.set_ylabel(r'$\frac{dR}{dt} (nm/t_{step})$',fontweight='bold',size=12)
cb=plt.colorbar(sc)
cb.set_label('Topological Class (n)',fontweight='bold',size=12)
# cb.set_label('Number Of Faces')
# plt.plot(xnew,poly1d_fn(xnew),linestyle='--',linewidth=2,c='black',label='y = %.3f x + %.3f' %(coef[0], coef[1]))
# plt.axhline(y=0, color='k', linestyle='-',linewidth=1)
# plt.axvline(x=0, color='k', linestyle='-',linewidth=1)
# plt.axvline(x=-coef[1]/coef[0], color='k', linestyle='-.',linewidth=1)
# ax.scatter(xa,np.array(G_tc))
# plt.axvline(x=-coef[1]/coef[0], color='k', linestyle='-.',linewidth=1)
# ax.legend(loc="best",prop={'size': 12})

# plt.ylim(-7500, 7500)
# sub_axes = plt.axes([.45, .5, .25, .25]) 
# sub_axes.scatter(xa,ya,c=view,cmap=cmap,norm=norm,s=50,marker='+')
# plt.ylim(-75, 75)
plt.tight_layout()
plt.show()
'''
dict = {'G_tc': G_tc, 'G_dv': G_dv, 'G_TL_l': G_TL_l,'G_TL_n':G_TL_n,'G_GB_n':G_GB_n,'G_mis':G_mis,'G_v':G_v,'G_a':G_a,'G_id':G_id,'G_tcf':G_tcf,'G_step':G_step,'G_dis':G_dis,'G_tmc':G_tmc,'G_af':G_af,'G_tcf_p':G_tcf_p,'G_a_mean':G_a_mean,'G_a_meanf':G_a_meanf,'G_GB_n_ext':G_GB_n_ext,'G_mis_w':G_mis_w}
df = pd.DataFrame(dict)
import os  
os.makedirs('E:/MicroTANA/G_data', exist_ok=True)  
df.to_csv('E:/MicroTANA/G_data/G_temp.csv')
'''