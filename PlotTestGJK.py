import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Arrow
from matplotlib.patches import Polygon
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import mpl_toolkits.mplot3d as a3
from scipy.spatial import ConvexHull
from GJK import *


xyz = np.array([
    [0,0,0],
    [1,0,0], 
    [1,1,0],
    [1,1,1]
    ], dtype=np.float) 

def minkowski_hull(poly_A, poly_B):
    na = poly_A.shape[0]
    nb = poly_B.shape[0]

    points = np.zeros([na * nb, poly_A.shape[1]])

    idx = 0
    for i in range(na):
        for j in range(nb):
            points[idx] = poly_A[i] - poly_B[j]
            idx += 1

    res = ConvexHull(points)
    vertices = points[res.vertices]
    return vertices

def set_limits(plt, polys):
    min = [np.min(poly, axis = 0) for poly in polys]
    max = [np.max(poly, axis = 0) for poly in polys]
    max = np.max(max, axis = 0)
    min = np.min(min, axis = 0)

    plt.xlim(min[0], max[0])
    plt.ylim(min[1], max[1])


def face_test():
    simplex = np.array(
            [[0.69228538, 1.1242913 ],
            [-0.75003475, -1.33943442],
            [-1.1899593 , -0.70945248]]
            )

    pc, best_simplex = face_closest_point(simplex)
    contains = contains_origin(simplex)


    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    patches = []
    polygon = Polygon(simplex, True)
    patches.append(polygon)

    circle = Circle([0,0], .01)
    patches.append(circle)

    collections = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

    colors = 100*np.random.rand(len(patches))
    collections.set_array(np.array(colors))

    ax.add_collection(collections)
    set_limits(plt, [simplex])

    plt.show()

    pass

#face_test()


np.random.seed(2039)

def plot_test_2d():

    total = 0
    while True:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')

        patches = []
    
        total = 1
        np.random.seed(total)
        poly_A = np.random.rand(3,2)*3 - 1
        poly_B = np.random.rand(3,2)*3 - 1

        if False:
            poly_A = np.array([
                [0,0],
                [2,1], 
                [0,2],
                ]) 

            poly_B = np.array([
                [3,0],
                [1,1], 
                [3,2],
                ]) 
            #poly_B -= 1
            #poly_A += 1


        print("Test {0}:".format(total))


        polygon = Polygon(poly_A, True)
        patches.append(polygon)

        polygon = Polygon(poly_B, True)
        patches.append(polygon)

        hull = minkowski_hull(poly_A, poly_B)
        polygon = Polygon(hull, True)
        patches.append(polygon)

        collections = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

        colors = 100*np.random.rand(len(patches))
        collections.set_array(np.array(colors))

        ax.add_collection(collections)

        print("Poly A")
        print(poly_A)
        print("Poly B")
        print(poly_B)
        print("Hull")
        print(hull)

        intersects = GJK(poly_A, poly_B)
        print(intersects)

        set_limits(plt, [hull, poly_A, poly_B])

        plt.title("{0} - {1}".format(total, intersects))

        plt.show()

plot_test_2d()
