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


# return convex hell, ordered count-clockwise
def convex_hull(points):
    res = ConvexHull(points)
    vertices = points[res.vertices]
    return vertices

def minkowski_hull(poly_A, poly_B):
    na = poly_A.shape[0]
    nb = poly_B.shape[0]

    points = np.zeros([na * nb, poly_A.shape[1]])

    idx = 0
    for i in range(na):
        for j in range(nb):
            points[idx] = poly_A[i] - poly_B[j]
            idx += 1

    return convex_hull(points)

def set_limits(plt, polys):
    min = [np.min(poly, axis = 0) for poly in polys]
    max = [np.max(poly, axis = 0) for poly in polys]
    max = np.max(max, axis = 0)
    min = np.min(min, axis = 0)

    plt.xlim(min[0], max[0])
    plt.ylim(min[1], max[1])


# return True if all points of polly are on the clockwise
# side of the line from b to a
def test_line_poly(point_a, point_b, poly):

    delta = point_b - point_a
    perp = np.array([delta[1], -delta[0]]) # rotate clockwise

    dp = perp.dot(point_a)    
    assert(perp.dot(point_b) >= dp - 1e-6)

    dists = poly.dot(perp)
    min = np.min(dists)

    return min >= dp

def test_poly_vs_segments(poly_A, poly_B):

    for i in range(poly_A.shape[0]):
        inext = (i + 1) % poly_A.shape[0] 

        if test_line_poly(poly_A[i], poly_A[inext], poly_B):
            return True

    return False

# simple test for convex poly intersection by check each segment a 
# a seperating axis. Returns true if there is an intersection:
def simple_2d_intersection_test(poly_A, poly_B):

    # test each segment of poly A to see it can be a separating axis:
    if test_poly_vs_segments(poly_A, poly_B):
        return False

    if test_poly_vs_segments(poly_B, poly_A):
        return False        

    return True


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
    
        total += 1
        np.random.seed(total)


        num_A = np.random.randint(3, 10)
        num_B = np.random.randint(3, 10)

        poly_A = convex_hull(np.random.rand(num_A,2)*2 - 1)
        poly_B = convex_hull(np.random.rand(num_B,2)*2 - 1)

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

        print("Poly A")
        print(poly_A)
        print("Poly B")
        print(poly_B)

        s_intersects = simple_2d_intersection_test(poly_A, poly_B)
        intersects = GJK(poly_A, poly_B)
        print(intersects)

        if intersects == s_intersects:
            continue
            pass

        fig, ax = plt.subplots()
        ax.set_aspect('equal')

        patches = []

        polygon = Polygon(poly_A, True)
        patches.append(polygon)

        polygon = Polygon(poly_B, True)
        patches.append(polygon)

        hull = minkowski_hull(poly_A, poly_B)
        print("Hull")
        print(hull)

        polygon = Polygon(hull, True)
        patches.append(polygon)

        collections = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

        colors = 100*np.random.rand(len(patches))
        collections.set_array(np.array(colors))

        ax.add_collection(collections)


        set_limits(plt, [hull, poly_A, poly_B])

        plt.title("{0} - {1}, {2}".format(total, intersects, s_intersects))

        plt.show()

plot_test_2d()
