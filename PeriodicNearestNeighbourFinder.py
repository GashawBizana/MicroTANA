
#SciPy-based kd-tree with periodic boundary conditions
#Copyright 2012 Patrick Varilly.
#Released under the scipy license

from periodic_kdtree import PeriodicCKDTree

# Find the index of the neighbours withn threshold distance in a periodic domain bounded by the edge length of the system
def PeriodicNeighbourIndex(bounds,data,ThrusholdDistance):
    Tree = PeriodicCKDTree(bounds, data)
    w2 = Tree.query_ball_point(data, ThrusholdDistance)
    return(w2)


