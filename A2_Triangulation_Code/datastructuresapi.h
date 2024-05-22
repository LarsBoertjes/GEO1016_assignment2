//
// Created by LarsB on 18/05/2024.
//

#ifndef A2_TRIANGULATION_CODE_DATASTRUCTURESAPI_H
#define A2_TRIANGULATION_CODE_DATASTRUCTURESAPI_H

/// Below are a few examples showing some useful data structures and APIs.

/// define a 2D vector/point
Vector2D b(1.1, 2.2);

/// define a 3D vector/point
Vector3D a(1.1, 2.2, 3.3);

/// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
Vector2D p = a.cartesian();

/// get the Homogeneous coordinates of p
Vector3D q = p.homogeneous();

/// define a 3 by 3 matrix (and all elements initialized to 0.0)
Matrix33 A;

/// define and initialize a 3 by 3 matrix
Matrix33 T(1.1, 2.2, 3.3,
           0, 2.2, 3.3,
           0, 0, 1);

/// define and initialize a 3 by 4 matrix
Matrix34 M(1.1, 2.2, 3.3, 0,
           0, 2.2, 3.3, 1,
           0, 0, 1, 1);

/// set first row by a vector
M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

/// set second column by a vector
M.set_column(1, Vector3D(5.5, 5.5, 5.5));

/// define a 15 by 9 matrix (and all elements initialized to 0.0)
Matrix W(15, 9, 0.0);
/// set the first row by a 9-dimensional vector
W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

/// get the number of rows.
int num_rows = W.rows();

/// get the number of columns.
int num_cols = W.cols();

/// get the the element at row 1 and column 2
double value = W(1, 2);

/// get the last column of a matrix
Vector last_column = W.get_column(W.cols() - 1);

/// define a 3 by 3 identity matrix
Matrix33 I = Matrix::identity(3, 3, 1.0);

/// matrix-vector product
Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

/// To access the value of an element.
double a = array[2];

/// define a 2D vector/point
Vector2D b(1.1, 2.2);

/// define a 3D vector/point
Vector3D c(1.1, 2.2, 3.3);

/// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
Vector2D p = c.cartesian();

/// get the Homogeneous coordinates of p
Vector3D q = p.homogeneous();

/// the length of a vector
double len = p.length();
/// the squared length of a vector
double sqr_len = p.length2();

/// the dot product of two vectors
double dot_prod = dot(p, q);

/// the cross product of two vectors
Vector cross_prod = cross(c, q);

/// normalize this vector
cross_prod.normalize();

#endif //A2_TRIANGULATION_CODE_DATASTRUCTURESAPI_H
