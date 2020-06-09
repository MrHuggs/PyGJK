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

# return convex hell, ordered count-clockwise
def convex_hull(points):
    res = ConvexHull(points)
    vertices = points[res.vertices]
    return vertices

# return the convex hull of hte Minkowski differnce.
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

def set_plot_limits(plt, polys):
    min = [np.min(poly, axis = 0) for poly in polys]
    max = [np.max(poly, axis = 0) for poly in polys]
    max = np.max(max, axis = 0)
    min = np.min(min, axis = 0)

    plt.xlim(min[0], max[0])
    plt.ylim(min[1], max[1])


# return True if all points of poly are on the clockwise
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

def seg_point_closest(a, b, p):

    delta = b - a
    p_delta = p - a

    dp = delta.dot(p_delta)

    if dp <= 0:
        return np.linalg.norm(p_delta)

    projection_2 = dp * dp / delta.dot(delta)
    p_delta_2 = p_delta.dot(p_delta)

    delta_2 = delta.dot(delta)

    if projection_2 >= delta_2:
        return np.linalg.norm(p - b)

    d2 = p_delta_2 - projection_2
    return np.math.sqrt(d2)


def poly_point_distance(poly, point):

    d = 1e+20
    for i in range(poly.shape[0]):
        inext = (i + 1) % poly.shape[0] 

        dc = seg_point_closest(poly[i], poly[inext], point)
        if (dc < d):
            d = dc

    return d

def poly_points_distance(poly_A, poly_B):
    d = 1e+20
    for b in poly_B:
        dc = poly_point_distance(poly_A, b)
        if dc < d:
            d = dc

    return d

def poly_poly_distance(poly_A, poly_B):
    dab = poly_points_distance(poly_A, poly_B)
    dba = poly_points_distance(poly_B, poly_A)

    if dab < dba:
        return dab
    return dba


def plot_test_2d():

    total = 0
    while True:
    
        total += 1
        np.random.seed(total)

        num_A = np.random.randint(3, 10)
        num_B = np.random.randint(3, 10)

        poly_A = convex_hull(np.random.rand(num_A,2)*2 - 1)
        poly_B = convex_hull(np.random.rand(num_B,2)*2)

        print("Test {0}:".format(total))

        print("Poly A")
        print(poly_A)
        print("Poly B")
        print(poly_B)

        s_intersects = simple_2d_intersection_test(poly_A, poly_B)
        intersects, gdist = GJK(poly_A, poly_B)
        print(intersects)

        if intersects == False:
            d = poly_poly_distance(poly_A, poly_B)
            print("Simple distance: {0}".format(d))
            print("GJK    distance: {0}".format(gdist))
            assert(np.allclose(d, gdist))
        

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


        set_plot_limits(plt, [hull, poly_A, poly_B])

        plt.title("{0} - {1}, {2}".format(total, intersects, s_intersects))

        plt.show()

plot_test_2d()
