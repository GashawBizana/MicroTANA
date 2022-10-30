import numpy as np
import trimesh
import pyvista as pv
import pymeshfix as mf

from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.pipeline import *
from ovito.vis import *
from ovito.qt_compat import QtCore

import open3d as o3d
import copy


def xyz_to_pcd(xyz):
    pcd=o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(xyz)
    
    return(pcd)

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
        
    return(pu_ac)

def Grain_mesh_ovito(GrainId,X,EdgeLength,grainId):
    
    data = DataCollection()
    cell = SimulationCell(pbc = (False, False, False))
    cell[:,0] = (EdgeLength[0],0,0)
    cell[:,1] = (0,EdgeLength[1],0)
    cell[:,2] = (0,0,EdgeLength[2])
    data.objects.append(cell)
    pos_u=[X[idx] for idx,i in enumerate(X) if GrainId[idx] ==grainId]
    pos=UnWrapCluster(xyz_to_pcd(np.array(pos_u)),EdgeLength)

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
        method = ConstructSurfaceModifier.Method.GaussianDensity, grid_resolution=50,
        radius_scaling = 1,
        isolevel = 0.9,only_selected=True)
    pipeline.modifiers.append(modifier_surface)
    # pipeline.modifiers.append(ConstructSurfaceModifier(radius = 7,only_selected=True))
    data1=pipeline.compute()
    mesh = data1.surfaces['surface']
    tri_mesh=trimesh.Trimesh(vertices=mesh.get_vertices(),faces=mesh.get_faces())
    return(tri_mesh)

def GrainBoundary_mesh_ovito(GrainBoundary_pos,EdgeLength):
    
    data = DataCollection()
    cell = SimulationCell(pbc = (False, False, False))
    cell[:,0] = (EdgeLength[0],0,0)
    cell[:,1] = (0,EdgeLength[1],0)
    cell[:,2] = (0,0,EdgeLength[2])
    data.objects.append(cell)
    pos_u=GrainBoundary_pos
    pos=UnWrapCluster(xyz_to_pcd(np.array(pos_u)),EdgeLength)

    gid=[0 for i in GrainBoundary_pos]
    particles = Particles()
    particles.create_property('Position', data=pos)
    particles.create_property('grainId', data=gid)
    data.objects.append(particles)

    # Create a new Pipeline with a StaticSource as data source:
    pipeline = Pipeline(source = StaticSource(data = data))

    modifier_select=ExpressionSelectionModifier(expression = "grainId==0")

    pipeline.modifiers.append(modifier_select)
    modifier_surface=ConstructSurfaceModifier(
        method = ConstructSurfaceModifier.Method.AlphaShape,
        radius = 3.5,smoothing_level=0,only_selected=True,transfer_properties=True)
    pipeline.modifiers.append(modifier_surface)
    # pipeline.modifiers.append(ConstructSurfaceModifier(radius = 7,only_selected=True))
    data1=pipeline.compute()
    mesh = data1.surfaces['surface']
    tri_mesh=trimesh.Trimesh(vertices=mesh.get_vertices(),faces=mesh.get_faces())
    return(tri_mesh)


def GrainProperty(InnerAtom,GrainId_list,X,EdgeLength):
    
    grain_property=[]
    
    for i in range (0,len(InnerAtom)):
        Grain_id=InnerAtom[i][0][0]
        mesh=Grain_mesh_ovito(GrainId_list,X,EdgeLength,Grain_id)
        grain_property.append([Grain_id,mesh])

    return(grain_property)

def GrainBoundaryProperty(GrainBoundary,EdgeLength):
    
    grain_boundary_property=[]
    
    for i in range (0,len(GrainBoundary)):
        Grain_boundary_id=GrainBoundary[i][0]
        Grain_boundary_pos=GrainBoundary[i][1:]
        mesh=GrainBoundary_mesh_ovito(Grain_boundary_pos,EdgeLength)
        grain_boundary_property.append([Grain_boundary_id,mesh])
    return(grain_boundary_property)


