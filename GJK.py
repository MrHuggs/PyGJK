import numpy as np
import itertools

def assert_is_array(a):
    assert(str(type(a)) == "<class 'numpy.ndarray'>")

def support(poly, direction):
    assert_is_array(poly)

    dot_products = poly.dot(direction)
    idx = np.argmax(dot_products)
    return idx

#
# find the closest point to the origin spanned by this simplex.
# return whether or not that point in obained in the simplex. 
#
def face_closest_point(simplex, tol = 1e-6):
    assert_is_array(simplex)

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

def closest_point(simplex, tol = 1e-6):
    assert_is_array(simplex)

    nump_points = len(simplex)

    best_p = None

    for i in range(1, nump_points + 1):
        for list in itertools.combinations(simplex, i):
            sub_simplex = np.vstack(list)
            in_face, p = face_closest_point(sub_simplex, tol)

            if in_face:
                d2 = p.dot(p)

                if best_p is None or d2 < best_d2:
                    best_d2 = d2
                    best_p = p
                    best_simplex = sub_simplex


    assert(not best_p is None)
    return best_p, best_simplex

def support_point(direction, poly_A, poly_B):
    assert_is_array(poly_A)
    assert_is_array(poly_B)

    idx_A = support(poly_A, direction)
    idx_B = support(poly_B, -direction)

    point = poly_A[idx_A] - poly_B[idx_B]
    return point


def GJK(poly_A, poly_B, epsilon = 1e-6):
    assert_is_array(poly_A)
    assert_is_array(poly_B)

    dimension = poly_A.shape[1]
    initial_dir = np.ones(dimension)
    initial_dir = np.array([0, -1])

    p = support_point(initial_dir, poly_A, poly_B)
    simplex = p.reshape([1, -1])

    while True:  

        next = support_point(-p, poly_A, poly_B)
        dp = p.dot(p) - next.dot(p)

        if dp < epsilon * epsilon:
            return False

        simplex = np.vstack([simplex, next])
        p, new_simplex = closest_point(simplex)

        pl2 = p.dot(p)

        if pl2 < epsilon * epsilon:
            return True

        simplex = new_simplex
        
    return True



