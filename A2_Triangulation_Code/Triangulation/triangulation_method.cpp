/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;
Matrix33 calculate_transformation_matrix(const std::vector<Vector2D> &points);
std::vector<Vector2D> normalize_points(const std::vector<Vector2D>& points, const Matrix33& T);
Matrix construct_matrix_W(const std::vector<Vector2D>& points_0_normalized, const std::vector<Vector2D>& points_1_normalized);
void printMatrix33(const Matrix33& M);
void printPoints(const std::vector<Vector2D> &points);
Matrix34 construct_matrix_M(const Matrix33& K, const Matrix33& R, const Vector3D& t);
Vector3D perform_linear_triangulation(const Vector2D &point1, const Vector2D &point2, const Matrix34& M, const Matrix34& M_prime);
int count_num_pos_z(const std::vector<Vector3D>& points3D);
void printPoints3D(const std::vector<Vector3D> &points);
Matrix34 concatenate(const Matrix33& R, const Vector3D& t);

/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        double s,                 /// input: the skew factor (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       triangulation(), or put them in one or multiple separate files.
    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).
    if (points_0.size() != points_1.size() || points_0.size() < 8) {
        std::cerr << "Error: mismatch in number of points or minimum required amount of 8 points is not met" << std::endl;
        return false;
    }


    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;

    // Compute transformation matrices
    Matrix33 T0 = calculate_transformation_matrix(points_0);
    Matrix33 T1 = calculate_transformation_matrix(points_1);

    // Normalize points
    std::vector<Vector2D> points_0_normalized = normalize_points(points_0, T0);
    std::vector<Vector2D> points_1_normalized = normalize_points(points_1, T1);

    printPoints(points_0_normalized);
    printPoints(points_1_normalized);

    // Construct matrix W
    Matrix W = construct_matrix_W(points_0_normalized, points_1_normalized);

    // Extract f matrix using SVD decomposition, W = USV.T, last column of V gives f
    Matrix U(W.rows(), W.rows());
    Matrix S(W.rows(), W.cols());
    Matrix V(W.cols(), W.cols());

    svd_decompose(W, U, S, V);

    Vector f = V.get_column(W.cols() - 1);

    // Convert f to F_estimate
    Matrix33 F_estimate(f[0], f[1], f[2],
                        f[3], f[4], f[5],
                        f[6], f[7], f[8]);

    printMatrix33(F_estimate);

    // Perform SVD decomposition again and change the smallest singular value to ensure rank 2
    Matrix U_est(F_estimate.rows(), F_estimate.rows());
    Matrix S_est(F_estimate.rows(), F_estimate.cols());
    Matrix V_est(F_estimate.cols(), F_estimate.cols());

    svd_decompose(F_estimate, U_est, S_est, V_est);

    // Ensure rank 2
    S_est(2,2) = 0;

    // Reconstruct fundamental matrix with rank 2
    Matrix33 Fnorm = U_est * S_est * V_est.transpose();

    // Denormalize fundamental matrix
    Matrix33 F = T1.transpose() * Fnorm * T0;

    // Step #2
    // construct K intrinsic matrix
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);

    // construct Essential matrix E
    Matrix E = K.transpose() * F * K;

    Matrix U_(E.rows(), E.rows());
    Matrix D_(E.rows(), E.cols());
    Matrix V_(E.cols(), E.cols());

    svd_decompose(E, U_, D_, V_);

    Matrix33 W_(0, -1, 0,
                1, 0, 0,
                0, 0, 1);

    Matrix33 R1 = determinant(U_ * W_ * V_.transpose()) * U_ * W_* V_.transpose();
    Matrix33 R2 = determinant(U_ * W_.transpose() * V_.transpose()) * U_ * W_.transpose() * V_.transpose();

    Vector3D t1 = U_.get_column(U_.cols() - 1);
    Vector3D t2 = -1 * U_.get_column(U_.cols() - 1);

    // Step 3
    // Construct M (all four options)
    Matrix33 R0 = Matrix33::identity();
    Vector3D t0 = (0.0, 0.0, 0.0);

    /*Matrix34 M = concatenate(K * R0, K * t0);
    Matrix34 M1_prime = concatenate(K * R1, K * t1);
    Matrix34 M2_prime = concatenate(K * R2, K * t1);
    Matrix34 M3_prime = concatenate(K * R1, K * t2);
    Matrix34 M4_prime = concatenate(K * R2, K * t2);

    Matrix34 M = construct_matrix_M(K, R0, t0);

    Matrix34 M1_prime = construct_matrix_M(K, R1, t1);
    Matrix34 M2_prime = construct_matrix_M(K, R2, t1);
    Matrix34 M3_prime = construct_matrix_M(K, R1, t2);
    Matrix34 M4_prime = construct_matrix_M(K, R2, t2);*/

    Matrix34 M = concatenate(K * R0, K * t0);

    // Get 4 M options for p1
    Matrix34 M1 = concatenate(K * R1, K * t1);
    Matrix34 M2 = concatenate(K * R1, K * t2);
    Matrix34 M3 = concatenate(K * R2, K * t1);
    Matrix34 M4 = concatenate(K * R2, K * t2);

    // linear triangulation method
    std::vector<Vector3D> points11, points12, points21, points22;

    for (int i = 0; i < points_0.size(); ++i) {
        points11.push_back(perform_linear_triangulation(points_0[i], points_1[i], M, M1));
        points12.push_back(perform_linear_triangulation(points_0[i], points_1[i], M, M2));
        points21.push_back(perform_linear_triangulation(points_0[i], points_1[i], M, M3));
        points22.push_back(perform_linear_triangulation(points_0[i], points_1[i], M, M4));
    }

    int count11 = count_num_pos_z(points11);
    int count12 = count_num_pos_z(points12);
    int count21 = count_num_pos_z(points21);
    int count22 = count_num_pos_z(points22);

    // Set t and R to the highest positive Z count
    int max_count = std::max({count11, count12, count21, count22});
    if (max_count == count11) {
        points_3d = points11;
        R = R1;
        t = t1;
    } else if (max_count == count12) {
        points_3d = points12;
        R = R1;
        t = t2;
    } else if (max_count == count21) {
        points_3d = points21;
        R = R2;
        t = t1;
    } else if (max_count == count22) {
        points_3d = points22;
        R = R2;
        t = t2;
    }

    std::cout << "Correct Matrix R" << std::endl;
    printMatrix33(R);

    std::cout << t << std::endl;

    return true;
}

Matrix33 calculate_transformation_matrix(const std::vector<Vector2D> &points) {
    size_t N = points.size();
    if (N == 0) return Matrix33(); // Return an identity matrix if no points

    // Calculate centroid
    double sumX = 0, sumY = 0;
    for (const auto& p : points) {
        sumX += p.x();
        sumY += p.y();
    }
    double centroidX = sumX / N;
    double centroidY = sumY / N;

    // Calculate scale factor
    double sumDistSquared = 0;
    for (const auto& p : points) {
        double dx = p.x() - centroidX;
        double dy = p.y() - centroidY;
        sumDistSquared += std::sqrt(dx * dx + dy * dy);
    }
    double avgDist = sumDistSquared / N;
    double scale = std::sqrt(2) / avgDist;

    // Create transformation matrix
    Matrix33 T;
    T(0, 0) = scale;
    T(0, 1) = 0;
    T(0, 2) = -scale * centroidX;
    T(1, 0) = 0;
    T(1, 1) = scale;
    T(1, 2) = -scale * centroidY;
    T(2, 0) = 0;
    T(2, 1) = 0;
    T(2, 2) = 1;

    return T;
}


std::vector<Vector2D> normalize_points(const std::vector<Vector2D>& points, const Matrix33& T) {
    std::vector<Vector2D> normalized_points;
    normalized_points.reserve(points.size());

    for (const Vector2D &p : points) {
        Vector3D p_homogeneous = p.homogeneous();
        Vector3D transformed = T * p_homogeneous;
        Vector2D normalized_p = transformed.cartesian();
        normalized_points.push_back(normalized_p);
    }

    return normalized_points;
}


Matrix construct_matrix_W(const std::vector<Vector2D>& points_0_normalized, const std::vector<Vector2D>& points_1_normalized) {
    Matrix W(points_0_normalized.size(), 9, 0.0);

    for (size_t i = 0; i < points_0_normalized.size(); ++i) {
        const auto& p0 = points_0_normalized[i];
        const auto& p1 = points_1_normalized[i];

        std::vector<double> row = {
                p1.x() * p0.x(),
                p1.x() * p0.y(),
                p1.x(),
                p1.y() * p0.x(),
                p1.y() * p0.y(),
                p1.y(),
                p0.x(),
                p0.y(),
                1.0
        };
        W.set_row(i, row);
    }

    return W;
}

void printMatrix33(const Matrix33& M) {
    for (int i = 0; i < M.rows(); ++i) {
        for (int j = 0; j < M.cols(); ++j) {
            std::cout << M(i, j) << " ";
        }
        std::cout << std::endl;
    }
}


void printPoints(const std::vector<Vector2D> &points) {
    for (const auto &point : points) {
        std::cout << "Point: (" << point.x() << ", " << point.y() << ")" << std::endl;
    }
}

void printPoints3D(const std::vector<Vector3D> &points) {
    for (const auto &point : points) {
        std::cout << "Point: (" << point.x() << ", " << point.y() << ", " << point.z() << ")" << std::endl;
    }
}

Matrix34 construct_matrix_M(const Matrix33& K, const Matrix33& R, const Vector3D& t) {
    Matrix34 M;

    // Compute K * R
    Matrix33 KR = K * R;

    // Set the first three columns of M to K * R
    M.set_column(0, Vector3D(KR(0,0), KR(1,0), KR(2,0)));
    M.set_column(1, Vector3D(KR(0,1), KR(1,1), KR(2, 1)));
    M.set_column(2, Vector3D(KR(0,2), KR(1,2), KR(2,2)));

    // Compute K * t
    Vector3D Kt = K * t;

    // Set the fourth column of M to K * t
    M.set_column(3, Kt);

    return M;
}


Matrix34 concatenate(const Matrix33& R, const Vector3D& t) {
    Matrix34 P;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            P[i][j] = R[i][j];
        }
        P[i][3] = t[i];
    }
    return P;
}

Vector3D perform_linear_triangulation(const Vector2D &point1, const Vector2D &point2, const Matrix34& M, const Matrix34& M_prime) {
    double x = point1[0];
    double y = point1[1];

    double x_prime = point2[0];
    double y_prime = point2[1];

    Vector4D m1 = M.get_row(0);
    Vector4D m2 = M.get_row(1);
    Vector4D m3 = M.get_row(2);

    Vector4D m1_prime = M_prime.get_row(0);
    Vector4D m2_prime = M_prime.get_row(1);
    Vector4D m3_prime = M_prime.get_row(2);

    // Construct Matrix A
    Matrix A(4, 4);
    A(0, 0) = x * m3[0] - m1[0];
    A(0, 1) = x * m3[1] - m1[1];
    A(0, 2) = x * m3[2] - m1[2];
    A(0, 3) = x * m3[3] - m1[3];

    A(1, 0) = y * m3[0] - m2[0];
    A(1, 1) = y * m3[1] - m2[1];
    A(1, 2) = y * m3[2] - m2[2];
    A(1, 3) = y * m3[3] - m2[3];

    A(2, 0) = x_prime * m3_prime[0] - m1_prime[0];
    A(2, 1) = x_prime * m3_prime[1] - m1_prime[1];
    A(2, 2) = x_prime * m3_prime[2] - m1_prime[2];
    A(2, 3) = x_prime * m3_prime[3] - m1_prime[3];

    A(3, 0) = y_prime * m3_prime[0] - m2_prime[0];
    A(3, 1) = y_prime * m3_prime[1] - m2_prime[1];
    A(3, 2) = y_prime * m3_prime[2] - m2_prime[2];
    A(3, 3) = y_prime * m3_prime[3] - m2_prime[3];

    Matrix U(A.rows(), A.rows());
    Matrix D(A.rows(), A.cols());
    Matrix V(A.cols(), A.cols());

    svd_decompose(A, U, D, V);

    Vector4D P_homogeneous = V.get_column(V.cols() - 1);
    Vector3D P(P_homogeneous[0] / P_homogeneous[3], P_homogeneous[1] / P_homogeneous[3], P_homogeneous[2] / P_homogeneous[3]);

    return P;
}

int count_num_pos_z(const std::vector<Vector3D>& points3D) {
    int count = 0;

    for (const Vector3D& pt : points3D) {
        if (pt.z() > 0) {
            count++;
        }
    }

    return count;
}