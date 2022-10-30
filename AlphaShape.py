from scipy.spatial import Delaunay
import numpy as np
import math



def alpha_shape(points, alpha, only_outer=True):
    """
    Compute the alpha shape (concave hull) of a set of points.
    :param points: np.array of shape (n, 3) points.
    :param alpha: alpha value.
    :param only_outer: boolean value to specify if we keep only the outer border
    or also inner edges.
    :return: set of (i,j) pairs representing edges of the alpha-shape. (i,j) are
    the indices in the points array.
    """
    assert np.shape(points)[0] > 3, "Need at least four points"

    def add_edge(edges, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the set already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            if only_outer:
                # if both neighboring triangles are in shape, it's not a boundary edge
                if (j, i) in edges:
                    edges.remove((j, i))
            return
        edges.add((i, j))

    tri = Delaunay(points)
    edges = set()
    # Loop over triangles:
    # ia, ib, ic, id = indices of corner points of the tetrahedron
    print(tri.vertices.shape)
    for ia, ib, ic, id in tri.vertices:
        pa = points[ia]
        pb = points[ib]
        pc = points[ic]
        pd = points[id]

        # Computing radius of tetrahedron Circumsphere
        # http://mathworld.wolfram.com/Circumsphere.html

        pa2 = np.dot(pa, pa)
        pb2 = np.dot(pb, pb)
        pc2 = np.dot(pc, pc)
        pd2 = np.dot(pd, pd)

        a = np.linalg.det(np.array([np.append(pa, 1), np.append(pb, 1), np.append(pc, 1), np.append(pd, 1)]))

        Dx = np.linalg.det(np.array([np.array([pa2, pa[1], pa[2], 1]),
                                     np.array([pb2, pb[1], pb[2], 1]),
                                     np.array([pc2, pc[1], pc[2], 1]),
                                     np.array([pd2, pd[1], pd[2], 1])]))

        Dy = - np.linalg.det(np.array([np.array([pa2, pa[0], pa[2], 1]),
                                       np.array([pb2, pb[0], pb[2], 1]),
                                       np.array([pc2, pc[0], pc[2], 1]),
                                       np.array([pd2, pd[0], pd[2], 1])]))

        Dz = np.linalg.det(np.array([np.array([pa2, pa[0], pa[1], 1]),
                                     np.array([pb2, pb[0], pb[1], 1]),
                                     np.array([pc2, pc[0], pc[1], 1]),
                                     np.array([pd2, pd[0], pd[1], 1])]))

        c = np.linalg.det(np.array([np.array([pa2, pa[0], pa[1], pa[2]]),
                                    np.array([pb2, pb[0], pb[1], pb[2]]),
                                    np.array([pc2, pc[0], pc[1], pc[2]]),
                                    np.array([pd2, pd[0], pd[1], pd[2]])]))

        circum_r = math.sqrt(math.pow(Dx, 2) + math.pow(Dy, 2) + math.pow(Dz, 2) - 4 * a * c) / (2 * abs(a))
        if circum_r < alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, id)
            add_edge(edges, id, ia)
            add_edge(edges, ia, ic)
            add_edge(edges, ib, id)
    return edges


# Plotting the output
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# creating figure
fig = plt.figure()
ax = Axes3D(fig)
  
for i in range(0,len(E)):
    p1=xyz[E[i][0]]
    p2=xyz[E[i][1]]
    x1=p1[0]
    x2=p2[0]
    y1=p1[1]
    y2=p2[1]
    z1=p1[2]
    z2=p2[2]
    ax.plot(p1,p2, color='green')
  
# setting title and labels
ax.set_title("3D plot")
ax.set_xlabel('x-axis')
ax.set_ylabel('y-axis')
ax.set_zlabel('z-axis')
  
# displaying the plot
plt.show()


  
# creating random dataset
xs = [14, 24, 43, 47, 54, 66, 74, 89, 12,
      44, 1, 2, 3, 4, 5, 9, 8, 7, 6, 5]
  
ys = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 6, 3,
      5, 2, 4, 1, 8, 7, 0, 5]
  
zs = [9, 6, 3, 5, 2, 4, 1, 8, 7, 0, 1, 2, 
      3, 4, 5, 6, 7, 8, 9, 0]
  

# creating the plot
ax.scatter(xs, ys, zs, color='green')
  
# setting title and labels
ax.set_title("3D plot")
ax.set_xlabel('x-axis')
ax.set_ylabel('y-axis')
ax.set_zlabel('z-axis')
  
# displaying the plot
plt.show()