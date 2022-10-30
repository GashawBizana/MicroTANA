import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
c1=1.51
Nc=13.3973
n=np.arange(5, 39, 1)
n_n=6-12/n
x_n=2*np.arctan((4*(np.sin(np.pi/n_n))**2-1)**(1/2))

G_h_1=((3*((((n-2)*np.tan(np.pi/n_n)))**(2/3))*(np.tan(x_n/2))**(1/3))/(2**(1/3)))*(np.pi/3-x_n)

G_h_2=(6/(2**(2/3)))*((3/(4*np.pi))**(1/3))*((np.tan(x_n/2))**(1/3))*((np.pi/3)-x_n)*(((n-2)*np.tan(np.pi/n_n))**(2/3))

# (3*((((n-2)*np.tan(np.pi/n_n)))**(2/3))*(np.tan(x_n/2))**(1/3))

# (np.pi/3-x_n)

dvdt_hill_1=2*G_h_1
dvdt_hill_2=2*G_h_2
#((3/(4*np.pi))**(1/3))*

# G_1=((3/(4*np.pi))**(1/3))*(np.pi/3-2*np.arctan(1.86*((n-1)**0.5)/(n-2)))

# G_2=((3/(4*np.pi))**(1/3))*( (5.35*(n**(2/3)))*((((n-2)/(2*((n-1)**0.5)))-((3/8)*G_1))**(-1/3)) )
# dvdt_mull=G_1*G_2

G_1=(np.pi/3-2*np.arctan(1.86*((n-1)**0.5)/(n-2)))

G_2=( (5.35*(n**(2/3)))*((((n-2)/(2*((n-1)**0.5)))-((3/8)*G_1))**(-1/3)) )

dvdt_mull=((3/(4*np.pi))**(1/3))*G_1*G_2

dvdt_rios=(3)*c1*(n**(1/2)-Nc**(1/2))
fig = plt.figure()
ax = fig.add_subplot()
ax.plot(n,dvdt_hill_1,'g.--',linewidth=1,c='red',label='Hilgenfeldt')
# ax.plot(n,dvdt_hill_2)
ax.plot(n,dvdt_mull,'b+--',linewidth=1,c='blue',label='Mullins')
ax.plot(n,dvdt_rios,'rx--',linewidth=1,c='green',label='Glicksman-Rios')

# plt.show()

F=np.array(G_GB_n)
ind_face=np.argsort(F)
F_ordered=F[ind_face]
G_dvOr=np.array(G_dv)*0.001*np.power(np.array(G_v)*0.001,-1/3)
G_dvOr_ordered=G_dvOr[ind_face]
srolv=-np.array(G_tcf)*0.1*(2*np.pi)*np.power(np.array(G_v)*0.001,-1/3)

srolv_app=-(np.array(G_tc)-np.array(G_TL_l)/6)*0.1*(2*np.pi)*np.power(np.array(G_v)*0.001,-1/3)

srolv_ordered=srolv[ind_face]
slope=np.array(G_dv)*0.001/((2*np.pi)*-np.array(G_tcf)*0.1)
slope_ordered=slope[ind_face]

srolv_app_ordered=srolv_app[ind_face]
mean_rate,num_f,mean_rate_srolv,mean_rate_srolv_app,srolv_std,srolv_app_std=[],[],[],[],[],[]
for i in range(min(F_ordered),max(F_ordered)+5):
    if i in (F_ordered):
        print(i)
        num_f.append(i)
        rate_i=[j for idx,j in enumerate(G_dvOr_ordered) if F_ordered[idx]==i and np.abs(G_dvOr_ordered[idx]) <1000]
        rate_i1=[j for idx,j in enumerate(srolv_ordered) if F_ordered[idx]==i and np.abs(G_dvOr_ordered[idx]) <1000]
        rate_i2=[j for idx,j in enumerate(srolv_app_ordered) if F_ordered[idx]==i and np.abs(G_dvOr_ordered[idx]) <1000]
        slop_i=[j for idx,j in enumerate(slope_ordered) if F_ordered[idx]==i and np.abs(G_dvOr_ordered[idx]) <1000]
        
        # mean_rate.append(np.mean(np.array(rate_i)/np.array(slop_i)))
        mean_rate.append(np.mean(np.array(rate_i)))
        mean_rate_srolv.append(np.mean(rate_i1))
        srolv_std.append(np.std(rate_i1))
        mean_rate_srolv_app.append(np.mean(rate_i2))
        srolv_app_std.append(np.std(rate_i2))

# ax.scatter(np.array(F_ordered),np.array(G_dvOr_ordered)/15,c='black', marker='+')
ax.scatter(np.array(F_ordered),srolv_ordered,c='black', marker='*',alpha=0.5,label='MachPerson-Srolovitz From Face')
ax.scatter(np.array(F_ordered),srolv_app_ordered,c='black', marker='+',alpha=0.5,label='MachPerson-Srolovitz assuming 120 deg dihedral angle')
# ax.plot(np.array(num_f),np.array(mean_rate)*(2*np.pi*1.51122595),marker='H',markersize= 5)
ax.plot(np.array(num_f),np.array(mean_rate_srolv),marker='s',markersize =5,c='orange',label='Mean of MachPerson-Srolovitz From Face')
# ax.errorbar(np.array(num_f),np.array(mean_rate_srolv),yerr=np.array(srolv_std),marker='s',markersize =5,c='orange')
ax.plot(np.array(num_f),np.array(mean_rate_srolv_app),marker='D',markersize =5,c='brown',label='Mean of MachPerson-Srolovitz assuming 120 deg dihedral angle')
# ax.errorbar(np.array(num_f),np.array(mean_rate_srolv_app),yerr=np.array(srolv_app_std),marker='D',markersize =5,c='brown')
# plt.ylim(-1000/15, 1000/15)
# plt.myButton.label.set_fontsize(12)
ax.set_xlabel(r'Topological Class (n)',fontweight='bold',size=12)
# ax.set_ylabel('dV/dt ($nm^3$/step)')
ax.set_ylabel(r'$\mathbf{\frac{dV}{dt} \frac{1}{V^{1/3}}  \frac{1}{m\gamma}}$',fontweight='bold',size=12)

ax.legend(loc="lower right",prop={'size': 6})
plt.show()
'''

def GrainsRecenter_plot(grain_mesh_m_1,EdgeLength):
    
    grain_mesh_n_0_cm=np.array([0,0,0])
    grain_mesh_m_1_cm=np.mean(grain_mesh_m_1.vertices,axis=0)
    
    diff_center=np.rint((grain_mesh_n_0_cm-grain_mesh_m_1_cm)/EdgeLength)    
    periodic_translation_factor=diff_center*EdgeLength
    grain_mesh_m_1.vertices= grain_mesh_m_1.vertices + periodic_translation_factor
    
    return(grain_mesh_m_1)

gb=[]
id=[0,2,1]
for k in [0,2]:
    EdgeLength=EdgeLength_for_all_timesteps[k]
    for i in range(0,len(grain_boundary_property_extended_for_all_timesteps[k])):
        if  grain_boundary_property_extended_for_all_timesteps[k][i][0][0] ==id[k] or id[k] == grain_boundary_property_extended_for_all_timesteps[k][i][0][0]  :
            gbi=copy.deepcopy(grain_boundary_property_extended_for_all_timesteps[k][i][1])
            gbi.remove_unreferenced_vertices()
            #trimesh.smoothing.filter_humphrey(gbi, alpha=1, beta=1, iterations=10, laplacian_operator=None)
            trimesh.smoothing.filter_taubin(gbi, lamb=0.5, nu=0.5, iterations=1, laplacian_operator=None)
            # gbi=GrainsRecenter_plot(gbi,EdgeLength)
            gbi.visual.vertex_colors = trimesh.visual.random_color()
            gb.append(gbi)        
        
combined = trimesh.util.concatenate(gb)
combined.show(viewer='gl') 
gm=gb[0]
gm.show(viewer='gl')
 
j=2
tlll=triple_line_length_for_all_timesteps[0]
gm=grain_property_for_all_timesteps[0][j][1]
gm=combined
fig2 = plt.figure(2)
ax3d = fig2.add_subplot(111, projection='3d')
for i in range(0,len(tlll)):
    if j in tlll[i][0]:
      TL=tlll[i][2]
      
      ax3d.plot(TL[:,0], TL[:,1], TL[:,2], 'o-')
if len(gm.triangles)>0:
    ax3d.plot_trisurf(gm.vertices[:, 0], gm.vertices[:,1], triangles=gm.faces, Z=gm.vertices[:,2],alpha=0.5) 
plt.show() 

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_trisurf(gm.vertices[:, 0], gm.vertices[:,1], triangles=gm.faces, Z=gm.vertices[:,2]) 
plt.show()

import numpy as np
import open3d as o3d
import trimesh
import copy
from scipy.spatial import distance



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
    return(pu_ac)

def GrainBoundaryRecenter(grain_mesh_n_0,gb,EdgeLength):
    grain_mesh_n_0_cm=np.mean(grain_mesh_n_0.vertices,axis=0)
    gbcm=np.mean(gb, axis=0)
    
    diff_center=np.rint((grain_mesh_n_0_cm-gbcm)/EdgeLength)    
    periodic_translation_factor=diff_center*EdgeLength
    gb= gb + periodic_translation_factor
    
    return(gb)

mesh=grain_property_for_all_timesteps[0][2][1]
EdgeLength=np.array([161.862, 162.025, 162.11 ])#EdgeLength_for_all_timesteps[0]
Xgb=np.array(GrainBoundary_for_all_timesteps[0][0][1:])

Xtl=[]

for i in range(0,len(TripleLine_for_all_timesteps[0])):
    if ((2 in TripleLine_for_all_timesteps[0][i][0]) and (1 in TripleLine_for_all_timesteps[0][i][0])):
        Xtl=Xtl+TripleLine_for_all_timesteps[0][i][1:]
Xqp=[]
for i in range(0,len(QuadraplePoint_for_all_timesteps[0])):
    if ((2 in QuadraplePoint_for_all_timesteps[0][i][0]) and (1 in QuadraplePoint_for_all_timesteps[0][i][0])):
        Xqp=Xqp+QuadraplePoint_for_all_timesteps[0][i][1:]
Xtl=np.array(Xtl)
Xqp=np.array(Xqp)

tol1=1*1.43
tol=1.5*1.43
m=copy.deepcopy(mesh)
m1=copy.deepcopy(mesh)

faces=m.faces
T=m.triangles
mask=[]
triangle=[]

pos_Xgb=GrainBoundaryRecenter(m,UnWrapCluster(xyz_to_pcd(Xgb),EdgeLength),EdgeLength)
pos_Xtl=GrainBoundaryRecenter(m,UnWrapCluster(xyz_to_pcd(Xtl),EdgeLength),EdgeLength)
pos_Xqp=GrainBoundaryRecenter(m,UnWrapCluster(xyz_to_pcd(Xqp),EdgeLength),EdgeLength)
triangle_centroid=np.mean(T,axis=1)
v_cdis=distance.cdist(triangle_centroid,pos_Xgb)
v_cdis1=distance.cdist(triangle_centroid,pos_Xtl)
v_cdis2=distance.cdist(triangle_centroid,pos_Xqp)
mask=np.logical_and(v_cdis.min(axis=1)<=tol,np.logical_and(v_cdis1.min(axis=1)>=tol1,v_cdis2.min(axis=1)>=tol1))
mask1=v_cdis.min(axis=1)<=tol
m.update_faces(mask)
m1.update_faces(mask1)
m.visual.vertex_colors = trimesh.visual.random_color()
m1.visual.vertex_colors = trimesh.visual.random_color()
mc = trimesh.util.concatenate([m,m1])
m1.show(viewer='gl')
m.show(viewer='gl')
print(m1.area)
print(m.area)

'''