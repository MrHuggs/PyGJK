import numpy as np
import itertools

def support(poly, direction):

    dot_products = poly.dot(direction)
    idx = np.argmax(dot_products)
    return idx

#
# find the closest point to the origin spanned by this simplex.
# return whether or not that point in obained in the simplex. 
#
def simplex_closest_point(simplex, tol = 1e-6):

    if simplex.shape[0] == 1:
        return True, simplex[0]

    x0 = simplex[0]
    d = simplex[1:] - x0

    D = d.dot(d.T)
    b = -d.dot(x0)

    try:
        a = np.linalg.solve(D, b)
    except np.linalg.LinAlgError:
        return False, x0

    min = np.min(a)
    sum = np.sum(a)

    in_face = min > -tol and sum < 1 + tol
    
    p = x0 + a.dot(d)
    return in_face, p

def simplex_check(simplex, tol = 1e-6):

    nump_points = len(simplex)

    best_p = None

    for i in range(1, nump_points + 1):
        for list in itertools.combinations(simplex, i):
            sub_simplex = np.vstack(list)
            in_face, p = simplex_closest_point(sub_simplex, tol)

            if in_face:
                d2 = p.dot(p)

                if best_p is None or d2 < best_d2:
                    best_d2 = d2
                    best_p = p
                    best_simplex = sub_simplex

    # Since we test 0-d simplices (points), there should alway be a best point:
    assert(not best_p is None)  
    return best_p, best_simplex

def support_point(direction, poly_A, poly_B):

    idx_A = support(poly_A, direction)
    idx_B = support(poly_B, -direction)

    point = poly_A[idx_A] - poly_B[idx_B]
    return point


# Test to convex shapes for intersection. Return true if they intersect
#
def GJK(poly_A, poly_B, epsilon = 1e-6):

    dimension = poly_A.shape[1]
    initial_dir = np.ones(dimension)

    p = support_point(initial_dir, poly_A, poly_B)
    simplex = p.reshape([1, -1])

    while True:  

        next = support_point(-p, poly_A, poly_B)
        dp = p.dot(p) - next.dot(p)

        if dp < epsilon * epsilon:
            return False, np.linalg.norm(p)

        simplex = np.vstack([simplex, next])
        p, new_simplex = simplex_check(simplex)

        pl2 = p.dot(p)

        if pl2 < epsilon * epsilon:
            return True, 0.

        simplex = new_simplex
        
    return True, 0.



