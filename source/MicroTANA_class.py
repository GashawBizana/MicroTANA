###### External Modules ######
import numpy as np
import copy



########## Import Local Modules ########
from Data_IO import  findAtomsGrain
from PeriodicNearestNeighbourFinder import PeriodicNeighbourIndex
from UnallocatedOrphanAtomGrainIdUpdater import UpdateGrainIdForUnallocatedOrphanAtoms
from UpdateX import update_x
from ModifiedGrainIdIdentifier import ModifiedGrainIdComputation
from GrainDataConstructor_modified import order_topologies
from Curvature_analysis_ovito import GrainMesh_and_Id,smoothing
from TripleLineLength_computation import TripleLineLength
from Microstructure_Topology_info import misorientation_quaternions,cubically_equivalent_quaternions,GrainBoundary_imc_sa_n_V,Grain_imc_sa_V


def GrainBoundary_mesh_modified(m_o,ll,id_grain0,id_grain1):
    m=copy.deepcopy(m_o)
    grainId=id_grain0
    Id_n=id_grain1
    faces=m.faces
    T=m.triangles
    mask=[]
    triangle=[]
    
    if True:
    
        for i in range(0,len(faces)):
            
        
            Test=[ (ll[faces[i][0]].count(Id_n)==1 and ll[faces[i][0]].count(grainId)==1 and ll[faces[i][0]].count(-1)==2) , (ll[faces[i][1]].count(Id_n)==1 and ll[faces[i][1]].count(grainId)==1 and ll[faces[i][1]].count(-1)==2), (ll[faces[i][2]].count(Id_n)==1 and ll[faces[i][2]].count(grainId)==1 and ll[faces[i][2]].count(-1)==2)]
            if Test.count(True) >=2:
                mask.append(True)
                triangle.append(faces[i])
            
            else :
                mask.append(False)
        m.update_faces(mask)  
    return(m,triangle)
        
class Grain:
    def __init__(self, grain_id, properties):
        self.grain_id = (grain_id,)
        self.properties = properties
    def compute_imc_surfacearea_volume_G(self):
        computed_properties=Grain_imc_sa_V(self.properties['mesh'])
        self.properties['surfacearea']=computed_properties['sa']
        self.properties['imc']=computed_properties['imc']
        self.properties['Volume']=computed_properties['V']
        
    def compute_topological_properties_G(self,M):   
        '''
        May need to have the information of GBs and TLs beforehand
        '''
        GB_linked_to_G=[GB_id for GB_id in M.grainboundary_id_list if set(self.grain_id).issubset(GB_id)]
        self.properties['n_GB']=len(GB_linked_to_G)
        self.properties['id_tl']=GB_linked_to_G
        
        
        tl_linked_to_G=[tl_id for tl_id in M.tripleline_id_list if set(self.grain_id).issubset(tl_id)]
        self.properties['n_tl']=len(tl_linked_to_G)
        self.properties['id_tl']=tl_linked_to_G
        qp_linked_to_G=[qp_id for qp_id in M.quadraplepoint_id_list if set(self.grain_id).issubset(qp_id)]
        self.properties['n_qp']=len(qp_linked_to_G)
        self.properties['id_qp']=qp_linked_to_G
        
        
        G_neighbors=[tuple([i for i in GB_id1 if i !=self.grain_id[0]]) for GB_id1 in GB_linked_to_G]
        
        self.properties['avg_n_GB_nG']=np.nanmean([len([GB_id for GB_id in M.grainboundary_id_list if set(gid_nn).issubset(GB_id)]) for gid_nn in G_neighbors])
        
        
        if set([ f'{tl_id}' for tl_id in tl_linked_to_G]).issubset(M.triplelines.keys()):
           
            self.properties['l_tl']=np.nansum([M.triplelines[f'{tl_id}'].properties['length'] for tl_id in tl_linked_to_G])  
        
        if set([ f'{GB_id}' for GB_id in GB_linked_to_G]).issubset(M.grainboundaries.keys()):
           
            self.properties['sa_from_GB']=np.nansum([M.grainboundaries[f'{GB_id}'].properties['surfacearea'][GB_id.index(GB_id[0])] for GB_id in GB_linked_to_G])  
            
            self.properties['imc_from_GB']=np.nansum([M.grainboundaries[f'{GB_id}'].properties['imc'][GB_id.index(GB_id[0])] for GB_id in GB_linked_to_G])
            
    def compute_relative_rotation(self,M,grain_mapping_data, cubMis):
        
        self.properties['Realative_rotation']=cubMis[M.timestep][np.where(grain_mapping_data[M.timestep] == self.grain_id[0])[0][0]]

    
        
            
        
        
class GrainBoundary:
    def __init__(self, grainboundary_id, properties):
        self.grainboundary_id = grainboundary_id
        self.properties = properties
    def compute_misorientation(self,misorientaion_data,grain_mapping_data):
        q0,q1,q2,q3=misorientaion_data
        q1_i,q2_i=misorientation_quaternions(self.properties['timestep'],self.grainboundary_id, q0,q1,q2,q3,grain_mapping_data)
        self.properties['misorient']=cubically_equivalent_quaternions(q1_i,q2_i)
            
    def compute_imc_surfacearea_volume_GB(self):
        computed_properties01=GrainBoundary_imc_sa_n_V(self.properties['mesh'][0])
        computed_properties10=GrainBoundary_imc_sa_n_V(self.properties['mesh'][1])
        self.properties['surfacearea']=(computed_properties01['sa'],computed_properties10['sa'])
        self.properties['imc']=(computed_properties01['imc'],computed_properties10['imc'])
        self.properties['normal']=(computed_properties01['GB_normal'],computed_properties10['GB_normal'])
        self.properties['Volume']=(computed_properties01['V'],computed_properties10['V']) #doesnt mean anything
        self.properties['imc_trimesh_positive']=(computed_properties01['imc_trimesh_positive'],computed_properties10['imc_trimesh_positive']) #doesnt mean anything
    def compute_topological_properties_GB(self,M):
        '''
        May need to have the information of TLs beforehand
        '''
        tl_linked_to_GB=[tl_id for tl_id in M.tripleline_id_list if set(self.grainboundary_id).issubset(tl_id)]
        self.properties['n_tl']=len(tl_linked_to_GB)
        self.properties['id_tl']=tl_linked_to_GB
        qp_linked_to_GB=[qp_id for qp_id in M.quadraplepoint_id_list if set(self.grainboundary_id).issubset(qp_id)]
        self.properties['n_qp']=len(qp_linked_to_GB)
        self.properties['id_qp']=qp_linked_to_GB
        
        if set([ f'{tl_id}' for tl_id in tl_linked_to_GB]).issubset(M.triplelines.keys()):
          self.properties['l_tl']=np.nansum([M.triplelines[f'{tl_id}'].properties['length'] for tl_id in tl_linked_to_GB])  
            

class TripleLine:
    def __init__(self, tripleline_id, properties):
        self.tripleline_id = tripleline_id
        self.properties = properties
        #def compute_misorientation(self,misorientaion_data):            
        
class Microstructure:
    def __init__(self,filename,timestep,ThresholdDistance=3.25,isPeriodic=True,atomic_radius=1.43):
        
        self.filename=filename
        self.timestep=timestep
        self.ThresholdDistance = ThresholdDistance
        self.isPeriodic = isPeriodic
        self.atomic_radius = atomic_radius
        self.grains={}
        self.grainboundaries={}
        self.triplelines={}

    def findAtomsGrain_c(self):
        self.X,self.GrainId,self.Image,self.H0,self.n_p,self.EdgeLength=findAtomsGrain(self.filename)
           

    def PeriodicNeighbourIndex_c(self):
       self.PeriodicNeighbourList_Index=PeriodicNeighbourIndex(self.EdgeLength,self.X,self.ThresholdDistance)
    
    def UpdateGrainIdForUnallocatedOrphanAtoms_c(self):
        UpdateGrainIdForUnallocatedOrphanAtoms(self.PeriodicNeighbourList_Index,self.GrainId,self.X,self.EdgeLength,self.isPeriodic)
    
    def ModifiedGrainIdComputation_c(self):
        self.ModifiedGrainIdListSorted=ModifiedGrainIdComputation(self.PeriodicNeighbourList_Index,self.GrainId,self.X,self.EdgeLength,self.isPeriodic)
        
    def Update_x_c(self):
        self.Xnew=update_x(self.X,self.ModifiedGrainIdListSorted,self.GrainId,self.EdgeLength,self.isPeriodic,self.PeriodicNeighbourList_Index)
    
    def Ordered_Topologies_c(self):
        self.Ordered_Topologies=order_topologies(self.ModifiedGrainIdListSorted,self.Xnew)
        self.n_grainboundaries=len(self.Ordered_Topologies['grain_boundary'])
        self.n_grains=len(self.Ordered_Topologies['inner_atom'])
        self.n_triplelines=len(self.Ordered_Topologies['triple_line'])
        self.n_quadraplepoints=len(self.Ordered_Topologies['quadraple_point'])
        self.grain_id_list=[i for i in range(self.n_grains)]
        self.grainboundary_id_list=[gb[0] for gb in self.Ordered_Topologies['grain_boundary']]
        self.tripleline_id_list=[tl[0] for tl in self.Ordered_Topologies['triple_line']]
        self.quadraplepoint_id_list=[qp[0] for qp in self.Ordered_Topologies['quadraple_point']]
    ############################################################

        
    def construct_grains(self,method=1,grainID_list=[]):
        if not grainID_list:
           ID_list=self.grain_id_list
        else:
            ID_list=grainID_list
        for id_i in ID_list:
            i=id_i[0]
            mid=[(self.ModifiedGrainIdListSorted[idx][2]).tolist() for idx,j in enumerate(self.Xnew) if self.GrainId[idx] ==i]
            gid=[self.GrainId[idx] for idx,j in enumerate(self.Xnew) if self.GrainId[idx] ==id_i]
            radius_list=[self.atomic_radius for idx,j in enumerate(self.Xnew) if self.GrainId[idx] ==id_i]
            InnerAtom_i=[self.Xnew[Idx].tolist() for Idx in range(len(self.Xnew)) if self.GrainId[Idx]==id_i]
            mesh_i,ll= GrainMesh_and_Id(InnerAtom_i,self.EdgeLength,mid,i,gid,radius_list,self.atomic_radius,method)
            #mesh_i.fix_normals()
            mesh_i=smoothing(copy.deepcopy(mesh_i),ll,0.01,100)
            grain_property={'timestep':self.timestep,'gid':i,'mesh':mesh_i,'v_property':ll}
            grain_i=Grain(id_i,grain_property)
            # add functionality to check if the grain already exists in the grain list
            self.grains[f'{id_i}']=grain_i
            
    def construct_grainboundaries(self,method=1,tolerance=[],grainboundaryID_list=[]):
        if not tolerance:
            
            tol=2.15*self.atomic_radius
            tol1=1.5*self.atomic_radius
        else:
            tol=tolerance[0] #maybe needed if another GB identification used 
            tol1=tolerance[1]
        if not grainboundaryID_list:
           ID_list=self.grainboundary_id_list
        else:
            ID_list=grainboundaryID_list
        for i in ID_list:
            id_grain0=(i[0],)
            id_grain1=(i[1],)
            mesh_grain0=self.grains[f'{id_grain0}'].properties['mesh']
            mesh_grain1=self.grains[f'{id_grain1}'].properties['mesh']
            ll_grain0=self.grains[f'{id_grain0}'].properties['v_property']
            ll_grain1=self.grains[f'{id_grain1}'].properties['v_property']
          
            # The below is better. Modify in the future
            #GrainBoundary_i=[gb[1] for gb in self.Ordered_Topologies['grain_boundary'] if gb[0]==i]
            grain_boundary_mesh_01,triangle01=GrainBoundary_mesh_modified(copy.deepcopy(mesh_grain0),ll_grain0,id_grain0,id_grain1)
            grain_boundary_mesh_10,triangle10=GrainBoundary_mesh_modified(copy.deepcopy(mesh_grain1),ll_grain1,id_grain1,id_grain0)
            
            grainboundary_property={'timestep':self.timestep,'gbid':i,'mesh':(grain_boundary_mesh_01,grain_boundary_mesh_10),'triangle_id':(triangle01,triangle10)}
            grainboundary_i=GrainBoundary(i,grainboundary_property)
            # add functionality to check if the grainboundary already exists in the grainboundary list
            self.grainboundaries[f'{i}']=grainboundary_i           
            
    
    def construct_triplelines(self,triplelineID_list=[]):
        
        if not triplelineID_list:
           ID_list=self.tripleline_id_list
        else:
            ID_list=triplelineID_list
        for i in ID_list:
            id_tripleline=i
            points_triplelines=np.array([tl[1] for tl in self.Ordered_Topologies['triple_line'] if tl[0]==id_tripleline][0])
            triple_line_length_i,interpolated_triple_line_path_i,start_end_i=TripleLineLength(points_triplelines,self.EdgeLength)
            tripleline_property={'timestep':self.timestep,'tlid':id_tripleline,'length':triple_line_length_i,'interpolated_path':interpolated_triple_line_path_i,'start_end':start_end_i}
            tripleline_i=TripleLine(id_tripleline,tripleline_property)
            self.triplelines[f'{id_tripleline}']=tripleline_i 