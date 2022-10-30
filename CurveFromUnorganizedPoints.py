import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.spatial as spatial
from sklearn.neighbors import KDTree
from scipy.spatial import Delaunay
import networkx as nx
import vg
import time
from scipy.optimize import curve_fit

def plot_tri_2(ax, points, tri):
    edges = collect_edges(tri)
    x = np.array([])
    y = np.array([])
    z = np.array([])
    for (i,j) in edges:
        x = np.append(x, [points[i, 0], points[j, 0], np.nan])      
        y = np.append(y, [points[i, 1], points[j, 1], np.nan])      
        z = np.append(z, [points[i, 2], points[j, 2], np.nan])
    ax.plot3D(x, y, z, color='g', lw='0.1')
    
    
def collect_edges(tri):
    edges = set()

    def sorted_tuple(a,b):
        return (a,b) if a < b else (b,a)
    # Add edges of tetrahedron (sorted so we don't add an edge twice, even if it comes in reverse order).
    for (i0, i1, i2, i3) in tri.simplices:
        edges.add(sorted_tuple(i0,i1))
        edges.add(sorted_tuple(i0,i2))
        edges.add(sorted_tuple(i0,i3))
        edges.add(sorted_tuple(i1,i2))
        edges.add(sorted_tuple(i1,i3))
        edges.add(sorted_tuple(i2,i3))
    return edges

    

def collect (Idx_P,H,A):
    
    A.append(Idx_P)
    W=[]
    W.append(1)
    #
    for node1, node2, data in mst_tree.edges(Idx_P, data=True):
        if data['weight'] <= H:
            A.append(int(node2))
            W.append((weight(data['weight'],H)))
    return(A,W)
            
def weight (r,H):
    
    w=np.exp(-np.square(r/H))
    
    return(w)

def f(xyz, A, B, C):
    """ Define Equation of a plan. """
    return A*xyz[:,0]+B*xyz[:,1]+C

def rms(z, zfit):
    return (np.sqrt(np.sum((z-zfit)**2)))


def planefit(xyz,w,guess=(10, 4, 2)):
    
    params, pcov = curve_fit(f, xyz[:,:2], xyz[:,2], guess,sigma=1/np.sqrt(np.array(w)), absolute_sigma=True)
    
    return params
def UnitNormal(parameter):
    coeff=np.array([parameter[0],parameter[1],1])
    n=coeff/np.linalg.norm(coeff)
    return n
    
    









# == RUN SCRIPT ========================================================================================================

start_time = time.perf_counter()

#Generate spiral points with noise 
total_rad = 10
z_factor = 3
noise = 0.06

num_sample_pts = 200
s_sample = np.linspace(0, total_rad, num_sample_pts)
x_sample = np.cos(s_sample) + noise * np.random.randn(num_sample_pts)
y_sample = np.sin(s_sample) + noise * np.random.randn(num_sample_pts)
z_sample = s_sample/z_factor + noise * np.random.randn(num_sample_pts)

# Make a [(x1, y1, z1), (x2, y2, z2)...] Numpy Array

points= np.vstack((x_sample, y_sample, z_sample )).T

#Compute Delaunay Triangulation

tri = Delaunay(points)
edge_lst = collect_edges(tri)
edge_lengths=[np.linalg.norm(points[e[0], :] - points[e[1], :]) for e in edge_lst]
G= nx.Graph((i, j, {'weight': dist}) for (i, j), dist in zip(edge_lst, edge_lengths))
mst_tree=nx.minimum_spanning_tree(G, weight='weight', algorithm='kruskal', ignore_nan=False)
mst_edges = nx.tree.minimum_spanning_edges(G, weight='weight',algorithm="kruskal", data=True,ignore_nan=False)
mst_edges_lst=list(mst_edges )

# Plotting
fig2 = plt.figure(2)
ax3d = fig2.add_subplot(111, projection='3d')

# Plot unordedered point cloud
ax3d.plot(points[:,0], points[:,1], points[:,2], 'm*')




































# xyz=np.array([[55.34482621, 13.81146467, 92.1027965 ],
#         [57.66349936, 13.75759136, 93.99737607],
#         [57.78505772, 15.11573352, 91.69995315],
#         [56.8020698 , 11.69741108, 88.88993841],
#         [55.07063198, 16.2096291 , 86.21431286],
#         [54.06368848, 13.37253895, 86.51713434],
#         [51.18918122, 15.23649075, 80.20538149],
#         [58.4611553 ,  9.65132697, 96.13593119],
#         [54.31344155, 15.95735618, 88.75165858],
#         [63.30439006, 13.1000615 , 95.21401162],
#         [48.85772098, 13.7850708 , 81.01398617],
#         [51.85265356, 15.02574483, 86.94364575],
#         [53.47402522, 15.11228238, 84.37614757],
#         [59.02411133, 11.4519432 , 93.97338379],
#         [56.53888219, 14.39841643, 96.2690235 ],
#         [60.90559522, 11.06810598, 95.99603026],
#         [58.74360449, 12.65165732, 96.38655325],
#         [58.88652863, 14.12908228, 98.6050286 ],
#         [59.14502225, 15.47241535, 96.25897268],
#         [61.08704252, 16.32223648, 98.29961336],
#         [61.37483316, 13.42101683, 97.54239755],
#         [56.42234155, 14.29225766, 87.90625493],
#         [60.6867578 , 13.5459219 , 94.6949354 ],
#         [59.75508013, 15.83682578, 93.60296244],
#         [52.18689859, 16.64207383, 82.50166964],
#         [61.62118712, 16.19099623, 95.57697591],
#         [58.65830321, 11.27034558, 98.66160499],
#         [55.99211235, 16.41734515, 83.50415788],
#         [54.67245146, 12.53064084, 83.83372751],
#         [57.56735333,  9.48379312, 92.91399494],
#         [57.59632663, 12.08573639, 91.5256849 ],
#         [49.35123821, 13.99056711, 78.35554428],
#         [51.28548911, 17.0916932 , 85.12704109],
#         [54.44859632, 13.33757395, 89.44662415]])