# PyGJK

A simple implementation of the GJK algorithm in Python.

While the GJK algorithm is elegant, when it is implemented for real-time use the result can be hard to understand.

Updating the simplex tends to be much larger than the core loop and has a lot of hairy 3d math. 

The update_simplex function in this version is more straightforward, but not something you would use in production. This version iterates over every subset of the vertices in the simplex, and does a traditional closest-point calculation to see if the origin is contained or what the closest point is.

numpy is used for linear algebra. A side-benefit is that it works in more than 3 dimensions.

# Usage

GJKExperiments.py

Test pyramids for intersection in arbitrary number of dimensions.

PlotTestGJK.py

Randomly creates 2d convex shapes and tests them for intersection. Draws the shapes and their Minkowski difference. If the shapes don’t overlap, draws the closest features in each shape.

# Commentary

There are a number of talks and tutorials on GJK, but I think they gloss over some important points:

Most of the work is in updating and checking the simplex. That is the function that will take the most time and code.

The termination condition, which is `|x|^2 – h_k(-x) < epsilon`, is not trivial and a key to making the algorithm work. In the original paper, it gets a separate proof in the appendix. The original paper was written in 1988, and the authors reference papers from 1976 and 1966 for the result.


# References

https://en.wikipedia.org/wiki/Gilbert%E2%80%93Johnson%E2%80%93Keerthi_distance_algorithm

The original paper: A fast procedure for computing the distance between complex objects in three-dimensional space, E.G. Gilbert ; D.W. Johnson; S.S. Keerthi 

Real-Time Collision Detection, Christer Ericson





