import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Arrow
from matplotlib.patches import Polygon
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import mpl_toolkits.mplot3d as a3
from GJK import *


xyz = np.array([
    [0,0,0],
    [1,0,0], 
    [1,1,0],
    [1,1,1]
    ], dtype=np.float) 


np.random.seed(2039)

def plot_test_2d():

    total = 0
    while True:
        fig, ax = plt.subplots()
        ax.set_aspect('equal')

        plt.xlim(-2, 2)
        plt.ylim(-2, 2)

        patches = []

        for i in range(5):
    
            np.random.seed(total)
            total += 1
            xyz = np.random.rand(3,2)*3 - 1

            print("Simplex {0}:".format(i))
            print(xyz)

            if contains_origin(xyz):
                print("Contains origin.")
                p = np.zeros(3)
            else:
                p = closest_point(xyz)
                print(p, np.linalg.norm(p))


            polygon = Polygon(xyz, True)
            patches.append(polygon)

            circle = Circle(p, .01)
            patches.append(circle)

            collections = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

            colors = 100*np.random.rand(len(patches))
            collections.set_array(np.array(colors))

            ax.add_collection(collections)

            ax.plot([0, p[0]],[0, p[1]])

        plt.show()

plot_test_2d()


def plot_test_3d():

    for i in range(1000):
    
        np.random.seed(i)
        xyz = np.random.rand(4,3)*3 - 1

        print("Simplex {0}:".format(i))
        print(xyz)

        if contains_origin(xyz):
            print("Contains origin.")
            p = np.zeros(3)
        else:
            p = closest_point(xyz)
            print(p, np.linalg.norm(p))

        indices = [0,1,2,0,1,3,0,2,3,1,2,3]

        axes = a3.Axes3D(pl.figure())

        tri = a3.art3d.Poly3DCollection(xyz[indices])
        tri.set_alpha(0.2)
        tri.set_color('grey')
        axes.add_collection3d(tri)

        axes.set_autoscale_on(True)
        axes.scatter(p[0],p[1],p[2])
        axes.scatter(0,0,0)

        max = np.max(xyz)
        min = np.min(xyz)
        axes.scatter(min, min, min)
        axes.scatter(max, max, max)

        pl.show()
