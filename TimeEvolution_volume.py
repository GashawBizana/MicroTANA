

grain_property=grain_property_for_all_timesteps[0]  
i=1
grain_i=grain_property[i][1]              

grain_i.show(viewer='gl')
grain_i.is_volume
print(grain_i.area)
print(grain_i.volume)
print(grain_i.integral_mean_curvature/np.pi)

'''
area = 25350

volume= 274625

mean width =195

'''



mesh=copy.deepcopy(cube_dn)

mesh_pv=pv.wrap(mesh)
meshfix = mf.MeshFix(mesh_pv)
meshfix.repair(verbose=False)
mesh_pv_repaired = meshfix.mesh
mesh=trimesh.Trimesh(mesh_pv_repaired.points,faces=mesh_pv_repaired.faces.reshape(-1, 4)[:, 1:])

face_adjacency_angles=mesh.face_adjacency_angles
edges_unique_length=mesh.edges_unique_length
face_adjacency_convex=mesh.face_adjacency_convex
mean_width=0
for i in range(0,len(edges_unique_length)):
    e_i=edges_unique_length[i]
    alpha_i=face_adjacency_angles[i]
    
    if face_adjacency_convex[i]==False:
        alpha_i=-alpha_i
    mean_width=mean_width+e_i*alpha_i
mean_width=(1/(2*np.pi))*mean_width
print(mean_width)


m=h0h1[0]
h0=m[0]
h1=m[1]
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
distnce=np.nanmean(h0n["distances"])

p = pv.Plotter()
p.add_mesh(h0n, scalars="distances", smooth_shading=True)
p.add_mesh(h1, color=True, opacity=0.75, smooth_shading=True)
p.show()

def GrainsRecenter(grain_mesh_n_0,grain_mesh_m_1,EdgeLength_step_i_plus_1):
    
    grain_mesh_n_0_cm=np.mean(grain_mesh_n_0.vertices,axis=0)
    grain_mesh_m_1_cm=np.mean(grain_mesh_m_1.vertices,axis=0)
    
    diff_center=np.rint((grain_mesh_n_0_cm-grain_mesh_m_1_cm)/EdgeLength_step_i_plus_1)    
    periodic_translation_factor=diff_center*EdgeLength_step_i_plus_1
    grain_mesh_m_1.vertices= grain_mesh_m_1.vertices + periodic_translation_factor
    
    return(grain_mesh_m_1)

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

El=EdgeLength#EdgeLength_for_all_timesteps[1]
m0=grain_property_for_all_timesteps[0][0][1]
h0=pv.wrap(GrainRepair(m0))
m1=grain_property_for_all_timesteps[1][5][1]
m1=GrainsRecenter(m0,m1,El)
h1=pv.wrap(GrainRepair(m1))

xyz=pv.wrap(np.array(InnerAtom_for_all_timesteps[0][0][1:]))
h0 = xyz.delaunay_3d(alpha=3.4, tol=0.001, offset=2.5, progress_bar=False)
h0= h0.extract_surface()
xyz1=pv.wrap(np.array(InnerAtom_for_all_timesteps[0][5][1:]))
h1 = xyz1.delaunay_3d(alpha=3.4, tol=0.001, offset=2.5, progress_bar=False)
h1= h1.extract_surface()
result = h0.boolean_intersection(h1)
pl = pv.Plotter()
_ = pl.add_mesh(h0, color='r', line_width=3)
_ = pl.add_mesh(h1, color='b', line_width=3) 
_ = pl.add_mesh(result, color='tan')
pl.show()


IA0=np.array(InnerAtom_for_all_timesteps[0][2][1:])

vtkpoints = PVGeo.points_to_poly_data(IA0)
bounds = vtkpoints.bounds
margin = 100
diameter = bounds[1] - bounds[0] + margin*2
n = 300
grid = pv.UniformGrid((n, n, n))
grid.origin = [bounds[0] - margin]*3
spacing = diameter/(n - 1)
grid.spacing = [spacing]*3

# interpolate the cloud's values onto the uniform grid
vox = grid.interpolate(vtkpoints, radius=spacing*2, progress_bar=True)

# extract the cell containing valid values
mask = vox['a'] > 0
vox_valid = vox.extract_points(mask, adjacent_cells=False)
vox_valid.plot()



xyz=np.array(InnerAtom_for_all_timesteps[0][0][1:])
xyz1=np.array(InnerAtom_for_all_timesteps[0][5][1:])

pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(xyz)

pcd1 = o3d.geometry.PointCloud()
pcd1.points = o3d.utility.Vector3dVector(xyz1)

mesh=o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd,7)
mesh.paint_uniform_color([0.5, 0.5, 0.5])
mesh1=o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd1,7)

o3d.visualization.draw_geometries([pcd,mesh], mesh_show_back_face=True)


'''


mid=[(ModifiedGrainIdListSorted[idx][2]).tolist() for idx,i in enumerate(X) if GrainId[idx] ==0]
for i in range(0,len(mid)):
    if len(mid[i])==1:
        mid[i]=[float(mid[i][0])]
        mid[i].append(-1.0)
        mid[i].append(-1.0)
        mid[i].append(-1.0)
    if len(mid[i])==2:
        mid[i]=[float(mid[i][0]),float(mid[i][1])]
        mid[i].append(-1.0)
        mid[i].append(-1.0)
    if len(mid[i])==3:
        mid[i]=[float(mid[i][0]),float(mid[i][1]),float(mid[i][2])]
        mid[i].append(-1.0)
    if len(mid[i])==4:
        mid[i]=[float(mid[i][0]),float(mid[i][1]),float(mid[i][2]),float(mid[i][3])]
        
data = DataCollection()
cell = SimulationCell(pbc = (False, False, False))
cell[:,0] = (EdgeLength[0],0,0)
cell[:,1] = (0,EdgeLength[1],0)
cell[:,2] = (0,0,EdgeLength[2])
data.objects.append(cell)
pos_u=[X[idx] for idx,i in enumerate(X) if GrainId[idx] ==0]
pos=UnWrapCluster(xyz_to_pcd(np.array(pos_u)),EdgeLength)

gid=[GrainId[idx] for idx,i in enumerate(X) if GrainId[idx] ==0]
particles = Particles()
particles.create_property('Position', data=pos)
particles.create_property('grainId', data=gid)
particles.create_property('la', data=mid)
data.objects.append(particles)
    
pipeline = Pipeline(source = StaticSource(data = data))

modifier_select=ExpressionSelectionModifier(expression = f"grainId=={0}")

pipeline.modifiers.append(modifier_select)
modifier_surface=ConstructSurfaceModifier(
    method = ConstructSurfaceModifier.Method.AlphaShape,
    radius = 4,smoothing_level=0,only_selected=True,transfer_properties=True)
pipeline.modifiers.append(modifier_surface)
# pipeline.modifiers.append(ConstructSurfaceModifier(radius = 7,only_selected=True))
data1=pipeline.compute()
mesh = data1.surfaces['surface']

l=mesh.vertices['la']
ll=[]
for i in l:
  ll.append(i.tolist())
  
m=trimesh.Trimesh(vertices=mesh.get_vertices(),faces=mesh.get_faces(), process=False)

mo3d=m.as_open3d
mo3dc=copy.deepcopy(mo3d)
triangle_vertices=np.asarray(mo3d.vertices)
verteices_to_remove=[]
for i in range(0,len(triangle_vertices)):
    if (1 not in ll[i].tolist() and 1 not in ll[i].tolist()):
        verteices_to_remove.append(i)
mo3dc.remove_unreferenced_vertices()

mo3dc.remove_vertices_by_index(verteices_to_remove)

mo3dc.paint_uniform_color([0.5, 0.5, 1])

o3d.visualization.draw_geometries([mo3dc,mo3d])
'''