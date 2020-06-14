###############################################################################
#
# Implementation of the Gilbert–Johnson–Keerthi distance algorithm.
# See:
#  * https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm
#  * the original paper: A fast procedure for computing the distance between complex objects 
#    in three-dimensional space,  E.G. Gilbert ; D.W. Johnson ; S.S. Keerthi 
#  * Real-Time Collision Detection, Christer Ericson
# 
#
import numpy as np
import itertools


###############################################################################
# Support function for a convex shape:
def support(poly, direction):
    dot_products = poly.dot(direction)
    idx = np.argmax(dot_products)
    return poly[idx]


###############################################################################
# Find the closest point to the origin in the affine hull of the simplex. 
# The affine hull is much larger than the convex hull because it allows negative coefficients,
# so the closest point might not be in the simplex. Return whether or not the closest
# point is in the simplex, and the closest point.
# 
def simplex_closest_point(simplex, tol=1e-6):
    if simplex.shape[0] == 1:
        return True, simplex[0]

    x0 = simplex[0]
    d = simplex[1:] - x0

    D = d.dot(d.T)
    b = -d.dot(x0)

    a = np.linalg.solve(D, b)

    min = np.min(a)
    sum = np.sum(a)

    in_face = min > -tol and sum < 1 + tol

    p = x0 + a.dot(d)
    return in_face, p


###############################################################################
# Given a list of points, find the point on or in the simplex closest to the origin (it will only
# be inside if the simplex contains the origin).
# If the closest point is not the origin, return closest point and the the smallest subset if the origin
# simplex so that the simplex made by the subset contains the closest point.
# The original points are assumed to be a in tuple, which is preserved (this is done so
# we can track the points in the original shapes).
#
def update_simplex(simplex_list, tol=1e-6):
    nump_points = len(simplex_list)

    best_p = None

    # To find the smallest sub-simplex with the closest point, we iterate over all possible
    # subsets, starting with the smallest, and keep the one that produces the smallest magnitude.
    #
    for i in range(1, nump_points + 1):
        for sub_tuple in itertools.combinations(simplex_list, i):
            sub_simplex = np.vstack([t[0] for t in sub_tuple])
            in_face, p = simplex_closest_point(sub_simplex, tol)

            if in_face:
                d2 = p.dot(p)

                if best_p is None or d2 < best_d2:
                    best_d2 = d2
                    best_p = p
                    best_sub_tuple = sub_tuple

    # Since we test 0-d simplices (points), there should alway be a best point:
    assert (best_p is not None)
    return best_p, list(best_sub_tuple)


def support_point(direction, poly_A, poly_B):
    a = support(poly_A, direction)
    b = support(poly_B, -direction)

    return a - b, (a, b)


###############################################################################
# Test to convex shapes for intersection. Return True if they intersect.
# If they don't intersect, return the closest distance between them, along
# with a list of tuple with points in each shape that form the closest features.
# Depending on the situation there might be duplications.
#
def GJK(poly_A, poly_B, epsilon=1e-6):
    dimension = poly_A.shape[1]

    initial_dir = np.ones(dimension)
    p, sups = support_point(initial_dir, poly_A, poly_B)
    simplex_list = [(p, sups)]

    while True:

        # This is the termination condition, which is not trivial to understand. See eqn 25
        # (and the proof in Appendix 1) in the original GJK paper for an
        # explanation:
        next, next_sups = support_point(-p, poly_A, poly_B)
        dp = p.dot(p) - next.dot(p)

        if dp < epsilon * epsilon:
            return False, np.linalg.norm(p), [t[1] for t in simplex_list]

        simplex_list.append((next, next_sups))
        p, new_simplex_list = update_simplex(simplex_list)

        pl2 = p.dot(p)

        if pl2 < epsilon * epsilon:
            return True, 0., None

        simplex_list = new_simplex_list
