
import numpy as np
from GJK import *



def face_closest_point_test():
    tests = [ 
        np.array([[1,0]]),
        np.array([[0,0.0000000001]]),
        np.array([[1,0], [0,1]]),
        np.array([[1,0,0], [0,1,0], [0,1,1]]),
        np.array([[1,0,0], [0,1,0], [0,0,1]]),
        np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0],[0,0,0,1]])
        ]
    cv = [
        [True, 1],
        [True, .0000000001],
        [True, 0.7071067811865476],
        [True, 0.7071067811865476],
        [True, 0.5773502691896257],
        [True, 0.5],
        ]

         
    for i, simp in enumerate(tests):
        in_face, p = simplex_closest_point(simp)
        print(simp)
        d = np.linalg.norm(p)
        print(in_face, p, d)
        ref_in_face, ref_d = cv[i]
        assert(ref_in_face == in_face)
        if in_face:
            assert(np.allclose(d, ref_d))

                       
face_closest_point_test()

def closest_point_test():
    tests = [ 
        np.array([[1,0,0,0], [0,1,0,0], [0,0,1,0],[0,0,0,1]])
        ]
    cv = [
        0.5773502691896257,
        ]

    for i, simp in enumerate(tests):
        p = closest_point(simp)

        print(simp)
        d = np.linalg.norm(p)
        print(p, d)
        ref_d = cv[i]
        assert(np.allclose(d, ref_d))

        pass

                       
closest_point_test()



def simplex_test():
    simp2 = np.array([[0,0], [1,0], [1,1]])

    assert(contains_origin(simp2 + 1) == False)
    assert(contains_origin(simp2 - .5))


    simp3 = np.array([[0,0,0], [1,0,0], [0,1,0],[0,0,1]])
    assert(contains_origin(simp2 + 1) == False)
    ac = simp2 - .5
    assert(contains_origin(ac))
    ac[[0,1]] = ac[[1,0]]
    assert(contains_origin(ac))

simplex_test()


def GJK_test():

    poly_A = np.array([
        [0,0],
        [1,0], 
        [1,1],
        [0,1]
        ]) 

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

    poly_B = poly_A + .5


    print(GJK(poly_A, poly_A + .5))
    print(GJK(poly_A, poly_A + 2))
