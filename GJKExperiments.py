###############################################################################
# Example of 3D intersection test
#
import sys
import numpy as np
from GJK import *


def make_cube(dim):
    n = np.power(2, dim)
    res = np.ones([n, dim])

    for i in range(n):
        for j in range(dim):
            if i & (1 << j) == 0:
                res[i, j] = -1

    return res


def make_pyramid(dim):
    t = make_cube(dim - 1)
    n = t.shape[0]

    new_col = np.zeros(n).reshape([-1, 1])
    t = np.hstack([t, new_col])

    new_row = np.zeros(dim)
    new_row[dim - 1] = 1

    res = np.vstack([new_row, t])
    return res


def GJK_test(dim):
    # 3D pyramid centered on the origin:
    shape_A = make_pyramid(dim)

    print(shape_A)

    offset = np.zeros(dim)
    offset[dim - 1] = 1
    shape_B = shape_A - offset * 2

    print(shape_B)

    # Test non-intersection
    print("Should not intersect:")
    print(GJK(shape_A, shape_B))

    # Test intersection
    shape_B = shape_A - offset * .5
    print("Should intersect:")
    print(GJK(shape_A, shape_B))


dim = 3
if len(sys.argv) > 1:
    dim = int(sys.argv[1])

print("Testing {0}-dimensional pyramids.".format(dim))

GJK_test(dim)
