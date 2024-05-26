
from periodic_kdtree import PeriodicCKDTree

#SciPy-based kd-tree with periodic boundary conditions
#Copyright 2012 Patrick Varilly.
#Released under the scipy license


def PeriodicNeighbourIndex(bounds,data,ThrusholdDistance):
    '''
    This function finds the index of the neighbours withn threshold distance in a periodic domain bounded by the edge length of the system.

    The function uses an updated verion of  SciPy-based kd-tree code that is written by Patrick Varilly and released under the scipy license. Please refer to https://github.com/patvarilly/periodic_kdtree/ for more details.
    
    '''
    Tree = PeriodicCKDTree(bounds, data)
    w2 = Tree.query_ball_point(data, ThrusholdDistance)
    return(w2)


