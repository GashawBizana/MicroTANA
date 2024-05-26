from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *

import numpy as np
from sklearn.cluster import DBSCAN
import trimesh
import copy
from scipy.spatial import distance
from scipy.sparse import coo_matrix
import random



def UnWrapCluster(pcd,EdgeLength):
    '''
    This function takes the position of a group of atoms in a periodic box
    and creat a cluster of this atom with a poition of their periodic image that gives the smallest non-periodic 
    distance. That is atoms are unrwaped so they are close to eachother
    
    '''

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


def GrainMesh_and_Id(pos_u,EdgeLength,mid,grainId,gid,radius_list,atomicradii,method):
    '''
    This function uses  "creat surface mesh " modifier of Ovito to creat a mesh represntation of grain boundaries of a grain and identify the mesh elements that corrospond to each grain boundary of a grain
    '''
    for i in range(0,len(mid)):
        if len(mid[i])==1:
            mid[i]=[float(mid[i][0])]
            mid[i].append(-1.0)
            mid[i].append(-1.0)
            mid[i].append(-1.0)
        if len(mid[i])==2:
            if mid[i][0]==gid[0]:
                mid[i]=[float(mid[i][0]),float(mid[i][1])]
            elif mid[i][1]==gid[0]:
                mid[i]=[float(mid[i][1]),float(mid[i][0])]
            mid[i].append(-1.0)
            mid[i].append(-1.0)
        if len(mid[i])==3:
            if mid[i][0]==gid[0]:
               mid[i]=[float(mid[i][0]),float(mid[i][1]),float(mid[i][2])]
            
            elif mid[i][1]==gid[0]:
                 mid[i]=[float(mid[i][1]),float(mid[i][0]),float(mid[i][2])]
            elif mid[i][2]==gid[0]:
                 mid[i]=[float(mid[i][2]),float(mid[i][0]),float(mid[i][1])]
            mid[i].append(-1.0)
        if len(mid[i])==4:
            if mid[i][0]==gid[0]:
               mid[i]=[float(mid[i][0]),float(mid[i][1]),float(mid[i][2]),float(mid[i][3])]
            elif mid[i][1]==gid[0]:
                 mid[i]=[float(mid[i][1]),float(mid[i][0]),float(mid[i][2]),float(mid[i][3])]
            elif mid[i][2]==gid[0]:
                  mid[i]=[float(mid[i][2]),float(mid[i][0]),float(mid[i][1]),float(mid[i][3])]
            elif mid[i][3]==gid[0]:
                  mid[i]=[float(mid[i][3]),float(mid[i][0]),float(mid[i][1]),float(mid[i][2])]  
                 
    data = DataCollection()
    cell = SimulationCell(pbc = (False, False, False))
    cell[:,0] = (EdgeLength[0],0,0)
    cell[:,1] = (0,EdgeLength[1],0)
    cell[:,2] = (0,0,EdgeLength[2])
    data.objects.append(cell)
    pos=UnWrapCluster(np.array(pos_u),EdgeLength)
    
    particles = Particles()
    particles.create_property('Radius', data=radius_list)
    particles.create_property('Position', data=pos)
    particles.create_property('grainId', data=gid)
    particles.create_property('la', data=mid)
    data.objects.append(particles)
    
        
    pipeline = Pipeline(source = StaticSource(data = data))
    
    modifier_select=ExpressionSelectionModifier(expression = f"grainId=={grainId}")
    
    pipeline.modifiers.append(modifier_select)
    # modifier_surface=ConstructSurfaceModifier(
    #     method = ConstructSurfaceModifier.Method.AlphaShape,
    #     radius =3.5,smoothing_level=0,only_selected=True,transfer_properties=True)
    
    if method==2 or method==3 or method==4 or method==5:
            modifier_surface=ConstructSurfaceModifier(
                method = ConstructSurfaceModifier.Method.GaussianDensity,
                grid_resolution=50,
                radius_scaling = 1,
                isolevel = atomicradii,only_selected=True,transfer_properties=True)
    elif method==1:
            modifier_surface=ConstructSurfaceModifier(
                method = ConstructSurfaceModifier.Method.AlphaShape,
                radius =3.25,smoothing_level=0,only_selected=True,transfer_properties=True)
    #4
    pipeline.modifiers.append(modifier_surface)
    # pipeline.modifiers.append(ConstructSurfaceModifier(radius = 7,only_selected=True))
    data1=pipeline.compute()
    mesh = data1.surfaces['surface']
    
    l=mesh.vertices['la']
    ll=[]
    for i in l:
      ll.append(i.tolist())
      
    m=trimesh.Trimesh(vertices=mesh.get_vertices(),faces=mesh.get_faces(), process=False)
    #m.fix_normals()
    
    return(m,ll)

def smoothing(m,ll,smoothfactor,iteration):
    '''
    This function contains several mesh smoothing method 
    '''
    
    def laplacian_calculation(m,equal_weight=False,pinned_vertices=[]):
        neighbors=m.vertex_neighbors
        
        # for i in pinned_vertices:
        #     neighbors[i]=[i]
        
        for idx,neghborlist in enumerate(neighbors):
            if idx in pinned_vertices:
                neighbors[idx]=[idx]
            else:
                neighbors[idx]=list(set(neghborlist)-set(pinned_vertices))
                if not neighbors[idx] :
                    neighbors[idx]=[idx]
                    
        vertices=m.vertices.view(np.ndarray)
        col=np.concatenate(neighbors)
        row=np.concatenate([[i]*len(n) for i,n in enumerate(neighbors)])
        
        if equal_weight:
            # equal weights for each neighbor
            data = np.concatenate([[1.0 / len(n)] * len(n)
                       for n in neighbors])
        else:
            # umbrella weights, distance-weighted
            # use dot product of ones to replace array.sum(axis=1)
            ones = np.ones(3)
            # the distance from verticesex to neighbors
            norms = [
                1.0 / np.maximum(1e-6, np.sqrt(np.dot(
                    (vertices[i] - vertices[n]) ** 2, ones)))
                for i, n in enumerate(neighbors)]
            # normalize group and stack into single array
            data = np.concatenate([i / i.sum() for i in norms])

        # create the sparse matrix
        matrix = coo_matrix((data, (row, col)),
                shape=[len(vertices)] * 2)
        
        return matrix
    
    def laplacian_calculation_tl(m,equal_weight=False,pinned_vertices=[]):
        neighbors=m.vertex_neighbors
        
        for idx,neghborlist in enumerate(neighbors):
            if idx in pinned_vertices:
                neighbors[idx]=[idx]
            else:
                neighbors[idx]=list(set(neghborlist)-set(pinned_vertices))
                if not neighbors[idx] :
                    neighbors[idx]=[idx]
                    

        vertices=m.vertices.view(np.ndarray)
        col=np.concatenate(neighbors)
        row=np.concatenate([[i]*len(n) for i,n in enumerate(neighbors)])
        
        if equal_weight:
            # equal weights for each neighbor
            data = np.concatenate([[1.0 / len(n)] * len(n)
                       for n in neighbors])
        else:
            # umbrella weights, distance-weighted
            # use dot product of ones to replace array.sum(axis=1)
            ones = np.ones(3)
            # the distance from verticesex to neighbors
            norms = [
                1.0 / np.maximum(1e-6, np.sqrt(np.dot(
                    (vertices[i] - vertices[n]) ** 2, ones)))
                for i, n in enumerate(neighbors)]
            # normalize group and stack into single array
            data = np.concatenate([i / i.sum() for i in norms])

        # create the sparse matrix
        matrix = coo_matrix((data, (row, col)),
                shape=[len(vertices)] * 2)
        
        return matrix
    
    
    
    TlOrQp=[i.count(-1)<2 for i in ll]
    NumberOfFaceVertices=len(TlOrQp)-sum(TlOrQp)
    
    NumberOfFaceVertices_10Percent=int(0.1*NumberOfFaceVertices)
    
    NumberOfTLOrQpVertices_90Percent=int(0.9*sum(TlOrQp))
    
    FaceVerticesIndex=[i for i in range(0,len(TlOrQp)) if TlOrQp[i]==False]
    
    for i in range(iteration):
        
        pinned_vertices=[idx for idx,i in enumerate(TlOrQp) if (i==False)]
        
        laplacian_matrix=laplacian_calculation_tl(m,equal_weight=False,pinned_vertices=pinned_vertices)
        trimesh.smoothing.filter_laplacian(m,
                      lamb=smoothfactor*0.5,
                      iterations=1,
                      implicit_time_integration=True,
                      volume_constraint=True,
                      laplacian_operator=laplacian_matrix)
    
        # FixedFaceVertices=random.sample(FaceVerticesIndex,k=NumberOfFaceVertices_10Percent)
        # FixedTLorQpVertices=random.sample([idx for idx,i in enumerate(TlOrQp) if (i==True)],k=NumberOfTLOrQpVertices_90Percent)
        pinned_vertices=[idx for idx,i in enumerate(TlOrQp) if (i==True)]
        # pinned_vertices=[idx for idx,i in enumerate(TlOrQp) if ((idx in FixedTLorQpVertices)or (idx in FixedFaceVertices))]
        # pinned_vertices=[idx for idx,i in enumerate(TlOrQp) if ((idx in FixedTLorQpVertices)or (idx in FixedFaceVertices))]
        
        laplacian_matrix=laplacian_calculation(m,equal_weight=False,pinned_vertices=pinned_vertices)
        trimesh.smoothing.filter_laplacian(m,
                      lamb=smoothfactor,
                      iterations=1,
                      implicit_time_integration=True,
                      volume_constraint=True,
                      laplacian_operator=laplacian_matrix)
        
    # for i in range(iteration):
        
    #     FixedFaceVertices=random.sample(FaceVerticesIndex,k=NumberOfFaceVertices_10Percent)
    #     FixedTLorQpVertices=random.sample([idx for idx,i in enumerate(TlOrQp) if (i==True)],k=NumberOfTLOrQpVertices_90Percent)
    #     # pinned_vertices=[idx for idx,i in enumerate(TlOrQp) if ((idx==True)or (idx in FixedFaceVertices))]
    #     pinned_vertices=[idx for idx,i in enumerate(TlOrQp) if ((idx in FixedTLorQpVertices)or (idx in FixedFaceVertices))]
        
    #     laplacian_matrix=laplacian_calculation(m,equal_weight=False,pinned_vertices=pinned_vertices)
    #     trimesh.smoothing.filter_laplacian(m,
    #                   lamb=smoothfactor,
    #                   iterations=1,
    #                   implicit_time_integration=False,
    #                   volume_constraint=True,
    #                   laplacian_operator=laplacian_matrix)
        
    return(m)
    
        
    

def test_for_clossnes(diff,tol,EdgeLength):
    Test0=[np.logical_or(np.logical_and(np.greater_equal(diff[:,0],-tol), np.less_equal(diff[:,0],tol)),np.logical_and(np.greater_equal(diff[:,0],EdgeLength[0]-tol),np.less_equal(diff[:,0],EdgeLength[0]+tol))),
           np.logical_or(np.logical_and(np.greater_equal(diff[:,1],-tol), np.less_equal(diff[:,1],tol)),np.logical_and(np.greater_equal(diff[:,1],EdgeLength[1]-tol),np.less_equal(diff[:,1],EdgeLength[1]+tol))),
           np.logical_or(np.logical_and(np.greater_equal(diff[:,2],-tol), np.less_equal(diff[:,2],tol)),np.logical_and(np.greater_equal(diff[:,2],EdgeLength[2]-tol),np.less_equal(diff[:,2],EdgeLength[2]+tol)))]
    

    return(np.any(np.all(Test0,axis=0)))

def GrainBoundaryRecenter(grain_mesh_n_0,gb,EdgeLength):
    grain_mesh_n_0_cm=np.mean(grain_mesh_n_0.vertices,axis=0)
    gbcm=np.mean(gb, axis=0)
    
    diff_center=np.rint((grain_mesh_n_0_cm-gbcm)/EdgeLength)    
    periodic_translation_factor=diff_center*EdgeLength
    gb= gb + periodic_translation_factor
    
    return(gb)

def GrainsRecenter(grain_mesh_n_0,grain_mesh_m_1,EdgeLength):
    
    grain_mesh_n_0_cm=np.mean(grain_mesh_n_0.vertices,axis=0)
    grain_mesh_m_1_cm=np.mean(grain_mesh_m_1.vertices,axis=0)
    
    diff_center=np.rint((grain_mesh_n_0_cm-grain_mesh_m_1_cm)/EdgeLength)    
    periodic_translation_factor=diff_center*EdgeLength
    grain_mesh_m_1.vertices= grain_mesh_m_1.vertices + periodic_translation_factor
    
    return(grain_mesh_m_1)



def GrainBoundary_mesh_ovito(m_o,ll,gbId,grainId,Xgb,EdgeLength,method,tol,tol1,mesh_j,Xtl):
    m=copy.deepcopy(m_o)
    mj=copy.deepcopy(mesh_j)
    
    Id_n=[i for i in gbId if i!=grainId][0]
    faces=m.faces
    T=m.triangles
    mask=[]
    triangle=[]
    
    if method==1:
    
        for i in range(0,len(faces)):
            
            # Test=[ (Id_n in ll[faces[i][0]]), (Id_n in ll[faces[i][1]]), (Id_n in ll[faces[i][2]])]
            # Test=[ (ll[faces[i][0]].count(Id_n)==1 and ll[faces[i][0]].count(-1)==2 ) , (ll[faces[i][1]].count(Id_n)==1 and ll[faces[i][0]].count(-1)==2), (ll[faces[i][2]].count(Id_n)==1 and ll[faces[i][0]].count(-1)==2)]
            Test=[ (ll[faces[i][0]].count(Id_n)==1 and ll[faces[i][0]].count(grainId)==1 and ll[faces[i][0]].count(-1)==2) , (ll[faces[i][1]].count(Id_n)==1 and ll[faces[i][1]].count(grainId)==1 and ll[faces[i][1]].count(-1)==2), (ll[faces[i][2]].count(Id_n)==1 and ll[faces[i][2]].count(grainId)==1 and ll[faces[i][2]].count(-1)==2)]
            if Test.count(True) >=2:
                mask.append(True)
                triangle.append(faces[i])
                
            # elif Test.count(True) ==1:
            #     Test1=[ ll[faces[i][0]].count(grainId)==1 and ll[faces[i][0]].count(-1)==3, ll[faces[i][1]].count(grainId)==1 and ll[faces[i][1]].count(-1)==3, ll[faces[i][2]].count(grainId)==1 and ll[faces[i][1]].count(-1)==3] 
                
            #     if Test1.count(True) >=1: # make this 2 if need to include triple line
            #         mask.append(True)
                    
            #     else:
            #         mask.append(False)
            
            else :
                mask.append(False)
        m.update_faces(mask)      
    elif method==2:
        
        for i in range(0,len(m.triangles)):
            v0,v1,v2=T[i][0],T[i][1],T[i][2]
            diff0,diff1,diff2=np.abs(v0-Xgb)%EdgeLength,np.abs(v1-Xgb)%EdgeLength,np.abs(v2-Xgb)%EdgeLength
            Test=[test_for_clossnes(diff0,tol,EdgeLength),test_for_clossnes(diff1,tol,EdgeLength),test_for_clossnes(diff2,tol,EdgeLength)]
            if Test.count(True) >=3:
                mask.append(True)
                triangle.append(faces[i])
                
            else :
                mask.append(False)
        m.update_faces(mask)
    elif method==3:
         v0,v1,v2=T[:,0],T[:,1],T[:,2]
         # test0,test1,test2=[test_for_clossnes(np.abs(vi-Xgb)%EdgeLength,tol,EdgeLength) for vi in v0],[test_for_clossnes(np.abs(vj-Xgb)%EdgeLength,tol,EdgeLength) for vj in v1],[test_for_clossnes(np.abs(vk-Xgb)%EdgeLength,tol,EdgeLength) for vk in v2]
         # Test=np.array([test0,test1,test2])
         # mask=np.all(Test,axis=0)
         mask=np.all([[test_for_clossnes(np.abs(vi-Xgb)%EdgeLength,tol,EdgeLength) for vi in T[:,0]],[test_for_clossnes(np.abs(vj-Xgb)%EdgeLength,tol,EdgeLength) for vj in T[:,1]],[test_for_clossnes(np.abs(vk-Xgb)%EdgeLength,tol,EdgeLength) for vk in T[:,2]]],axis=0)
         triangle=faces[np.where(mask==True)]
         m.update_faces(mask)
         
    elif method==4:
        triangle_centroid=np.mean(T,axis=1)
        if not (Xtl) ==False:
            pos_Xtl=GrainBoundaryRecenter(m,UnWrapCluster(Xtl,EdgeLength),EdgeLength)
            v_cdis_1=distance.cdist(triangle_centroid,pos_Xtl)
            pos_Xgb=GrainBoundaryRecenter(m,UnWrapCluster(Xgb,EdgeLength),EdgeLength)
            v_cdis=distance.cdist(triangle_centroid,pos_Xgb)
            mask=np.logical_and(v_cdis.min(axis=1)<=tol, v_cdis_1.min(axis=1)>=tol1)
        elif not (Xtl):
            pos_Xgb=GrainBoundaryRecenter(m,UnWrapCluster(Xgb,EdgeLength),EdgeLength)
            v_cdis=distance.cdist(triangle_centroid,pos_Xgb)
            mask=v_cdis.min(axis=1)<=tol
            
        triangle=faces[np.where(mask==True)]
        m.update_faces(mask)
        
        # for jj in range(0,2):
        #     if len(m.faces) >=10:
        #         broken = trimesh.repair.broken_faces(m)
        #         mask1=[(k not in broken) for k in range(len(m.faces))]
        #         m.update_faces(mask1)
    elif method==5:
        mj=GrainsRecenter(m,mj,EdgeLength)
        
        triangle_centroid=np.mean(T,axis=1)
        v_cdis=distance.cdist(triangle_centroid,np.mean(mj.triangles,axis=1))
        mask=v_cdis.min(axis=1)<=(tol)
        triangle=faces[np.where(mask==True)]
        m.update_faces(mask)
    
    return(m,triangle)


#def generate_grain_mesh(grain_id,inner_atom,EdgeLength,GrainID ):
    
    

def GrainBoundaryProperties_ovito(grain_boundary_list,grain_inner_atom_list,EdgeLength,X,GrainId,ModifiedGrainIdListSorted,TripleLine,QuadraplePoint,atomicradii,method):
    '''
    This function construct a list containing properties of a each grain, grain boundary, triple line and quadraple 
    junction in a given microstructure 
    '''
    number_of_grain_boundaries=len(grain_boundary_list)
    number_of_grains=len(grain_inner_atom_list)
    grain_boundary_property=[]
    grain_property=[]
    ll_list=[]
    tol=2.15*atomicradii
    tol1=1.5*atomicradii
    #1.25
    for i in range(0,number_of_grains):
        # InnerAtom_not_gb=grain_inner_atom_list[i][1:]
        mid=[(ModifiedGrainIdListSorted[idx][2]).tolist() for idx,j in enumerate(X) if GrainId[idx] ==i]
        gid=[GrainId[idx] for idx,j in enumerate(X) if GrainId[idx] ==i]
        radius_list=[atomicradii for idx,j in enumerate(X) if GrainId[idx] ==i]
        InnerAtom_i=[X[Idx].tolist() for Idx in range(len(X)) if GrainId[Idx]==i]
        
        mesh_i,ll= GrainMesh_and_Id(InnerAtom_i,EdgeLength,mid,i,gid,radius_list,atomicradii,method)
        mesh_i=smoothing(copy.deepcopy(mesh_i),ll,0.01,100)
        
        grain_property.append([i,mesh_i,ll])
        ll_list.append([i,ll])
    
    
    for i in range(0,number_of_grains):
        mesh_i=grain_property[i][1]
        ll=ll_list[i][1]
        
        for j in range (0, number_of_grain_boundaries):
            
            if i in grain_boundary_list[j][0]: #(i== grain_boundary_list[j][0][0] or i==grain_boundary_list[j][0][1]):
                Id_j=[idg for idg in grain_boundary_list[j][0] if idg!=i][0] # uncomment this for method 5
                GrainBoundary_i=[X[Idx].tolist() for Idx in range(len(X)) if (GrainId[Idx]==Id_j and i in (ModifiedGrainIdListSorted[Idx][2]).tolist() and len((ModifiedGrainIdListSorted[Idx][2]).tolist())==2)]
                
                TripleLine_i=[X[Idx].tolist() for Idx in range(len(X)) if (GrainId[Idx]==i and Id_j in (ModifiedGrainIdListSorted[Idx][2]).tolist() and len((ModifiedGrainIdListSorted[Idx][2]).tolist())>=3)]
                
                mesh_j=mesh_i#grain_property[Id_j][1]  Fix this if using method 5
                if method==1:
                    grain_boundary_mesh_ij,triangle=GrainBoundary_mesh_ovito(copy.deepcopy(mesh_i),ll,grain_boundary_list[j][0],i,np.array(grain_boundary_list[j][1:]),EdgeLength,method,tol,tol1,mesh_j,TripleLine_i)
                if method==4:
                    if not(GrainBoundary_i):
                        grain_boundary_mesh_ij,triangle=GrainBoundary_mesh_ovito(copy.deepcopy(mesh_i),ll,grain_boundary_list[j][0],i,np.array(grain_boundary_list[j][1:]),EdgeLength,method,tol,tol1,mesh_j,TripleLine_i)
                    else:
                        grain_boundary_mesh_ij,triangle=GrainBoundary_mesh_ovito(copy.deepcopy(mesh_i),ll,grain_boundary_list[j][0],i,np.array(GrainBoundary_i),EdgeLength,method,tol,tol1,mesh_j,TripleLine_i)
               
                
                grain_boundary_mesh_Trimesh= copy.deepcopy(grain_boundary_mesh_ij)
                

                grain_boundary_property.append([[i,[k for k in grain_boundary_list[j][0] if k!=i][0]],grain_boundary_mesh_Trimesh,triangle])
                
    return(grain_boundary_property,grain_property)



def GrainBoundaryPropertyExtended_ovito(GrainBoundary,grain_boundary_property_extended):
    '''
    This function construct a list containing properties of a each grain, grain boundary, triple line and quadraple 
    junction in a given microstructure 
    
    The only difference with GrainBoundaryProperties_ovito is that  this conatins GB mesh constructed by considering 
    each grain that share the GB
    '''
    grain_boundary_property=[]
    
    for i in range (0,len(GrainBoundary)):
        gb_Id=GrainBoundary[i][0]
        
        gb_meshes_extended=[j[1] for j in grain_boundary_property_extended if ((gb_Id[0] in j[0]) and (gb_Id[1] in j[0]))]
        
        grain_boundary_property.append([gb_Id,gb_meshes_extended])
        
    return(grain_boundary_property)
    
            
            
                
              
              
#### Code graveyard
'''    
        
gb=[]
for i in range(0,len(grain_boundary_property_for_all_timesteps[0])):
    if 70 in grain_boundary_property_for_all_timesteps[0][i][0]:
        gbi=grain_boundary_property_for_all_timesteps[0][i][1][0]
        gbi.visual.vertex_colors = trimesh.visual.random_color()
        gb.append(gbi)        
combined = trimesh.util.concatenate(gb)   
combined.show(viewer='gl') 
gb[0].show(viewer='gl')   

gb=[]
for i in range(0,len(grain_boundary_property_extended_for_all_timesteps[0])):
    if 3 == grain_boundary_property_extended_for_all_timesteps[0][i][0][0]:
        gbi=grain_boundary_property_extended_for_all_timesteps[0][i][1]
        gbi.visual.vertex_colors = trimesh.visual.random_color()
        gb.append(gbi)        
combined = trimesh.util.concatenate(gb)   
combined.show(viewer='gl') 
gb[0].show(viewer='gl')

fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot(111, projection='3d')
for i in range(0,len(xyz)):
    
    idx=i
    
    if (xyz[idx].shape)[0] >15:
        
        ax.plot(xyz[idx][:,0],xyz[idx][:,1],xyz[idx][:,2])

plt.show()

i=0
mid=[(ModifiedGrainIdListSorted[idx][2]).tolist() for idx,j in enumerate(X) if (i in (ModifiedGrainIdListSorted[idx][2]).tolist() and len((ModifiedGrainIdListSorted[idx][2]).tolist())>=1) ]


gid=[GrainId[idx] for idx,j in enumerate(X) if ( i in (ModifiedGrainIdListSorted[idx][2]).tolist() and len((ModifiedGrainIdListSorted[idx][2]).tolist())>=1)]



InnerAtom_i=[X[Idx].tolist() for Idx in range(len(X)) if ( i in (ModifiedGrainIdListSorted[Idx][2]).tolist() and len((ModifiedGrainIdListSorted[Idx][2]).tolist())>=1)]
mesh_i,ll= GrainMesh_and_Id(InnerAtom_i,EdgeLength,mid,i,gid)
mesh_i.show(viewer='gl')  

j=1
grain_boundary_mesh_ij=GrainBoundary_mesh_ovito(copy.deepcopy(mesh_i),ll,GrainBoundary[j][0],i)

grain_boundary_mesh_ij.show(viewer='gl')  


   

i=0
mid=[(ModifiedGrainIdListSorted[idx][2]).tolist() for idx,j in enumerate(X) if (i in (ModifiedGrainIdListSorted[idx][2]).tolist() and len((ModifiedGrainIdListSorted[idx][2]).tolist())>1) ]

gid=[GrainId[idx] for idx,j in enumerate(X) if ( i in (ModifiedGrainIdListSorted[idx][2]).tolist() and len((ModifiedGrainIdListSorted[idx][2]).tolist())>1)]

InnerAtom_i=[X[Idx].tolist() for Idx in range(len(X)) if ( i in (ModifiedGrainIdListSorted[Idx][2]).tolist() and len((ModifiedGrainIdListSorted[Idx][2]).tolist())>1)]
mesh_i,ll= GrainMesh_and_Id(InnerAtom_i,EdgeLength,mid,i,gid)
mesh_i.show(viewer='gl')  

j=1
grain_boundary_mesh_ij=GrainBoundary_mesh_ovito(copy.deepcopy(mesh_i),ll,GrainBoundary[j][0],i)

grain_boundary_mesh_ij.show(viewer='gl')  


i=0
mid=[(ModifiedGrainIdListSorted[idx][2]).tolist() for idx,j in enumerate(X) if (GrainId[idx] ==i and len((ModifiedGrainIdListSorted[idx][2]).tolist())>=1) ]
gid=[GrainId[idx] for idx,j in enumerate(X) if ( GrainId[idx] ==i and len((ModifiedGrainIdListSorted[idx][2]).tolist())>=1)]
InnerAtom_i=[X[Idx].tolist() for Idx in range(len(X)) if ( GrainId[Idx] ==i and len((ModifiedGrainIdListSorted[Idx][2]).tolist())>=1)]
mesh_i,ll= GrainMesh_and_Id(InnerAtom_i,EdgeLength,mid,i,gid)
mesh_i.show(viewer='gl') 

j=2
grain_boundary_mesh_ij=GrainBoundary_mesh_ovito(copy.deepcopy(mesh_i),ll,GrainBoundary[j][0],i)

grain_boundary_mesh_ij.show(viewer='gl') 

i=0
m=copy.deepcopy(grain_property_for_all_timesteps[0][i][1])
m.show(viewer='gl')
m1=copy.deepcopy(m)
idlist=copy.deepcopy(grain_property_for_all_timesteps[0][i][2])
m=smoothing(m,idlist,0.01,100)
m1.visual.vertex_colors = trimesh.visual.random_color()
combined = trimesh.util.concatenate([m1,m])
combined.show(viewer='gl')
m.show(viewer='gl')

import trimesh
import copy
def GrainsRecenter_plot(grain_mesh_m_1,EdgeLength):
    
    grain_mesh_n_0_cm=np.array([0,0,0])
    grain_mesh_m_1_cm=np.mean(grain_mesh_m_1.vertices,axis=0)
    
    diff_center=np.rint((grain_mesh_n_0_cm-grain_mesh_m_1_cm)/EdgeLength)    
    periodic_translation_factor=diff_center*EdgeLength
    grain_mesh_m_1.vertices= grain_mesh_m_1.vertices + periodic_translation_factor
    
    return(grain_mesh_m_1)

g=[]
k=0
EdgeLength=EdgeLength_for_all_timesteps[k]
ng=len(grain_property_for_all_timesteps[k])
for i in range(0,ng):
    m=copy.deepcopy(grain_property_for_all_timesteps[k][i][1])
    m=GrainsRecenter_plot(m,EdgeLength)
    m.visual.vertex_colors = trimesh.visual.random_color()
    g.append(m)
combined = trimesh.util.concatenate(g)
combined.show(viewer='gl')
    
'''
