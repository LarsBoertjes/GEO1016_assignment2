#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>
#include <cmath>
#include <numeric> // for std::accumulate

using namespace easy3d;

void print_vector(const Vector& p) {
    std::cout << "( ";
    for (int i = 0; i < p.size(); i++) {
        std::cout << p[i] << " ";
    }
    std::cout << ")" << std::endl;
}

void print_matrix(const Matrix& m) {
    int rows = m.rows();
    int cols = m.cols();
    std::cout << "---------------- Matrix ----------------" << std::endl;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << std::setw(10) << m.get(i, j) << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;
}

// Compute the distribution of 3D points in the two camera coordinate systems
void compute_point_distribution(
        const std::vector<Vector3D>& points3D,
        const Matrix34& M_0,
        const Matrix34& M_1
) {
    int front_of_camera_0 = 0;
    int front_of_camera_1 = 0;
    for (const Vector3D& P : points3D) {
        // Point in the first camera coordinate system
        Vector3D P_cam0 = M_0 * P.homogeneous();
        if (P_cam0.z() > 0) {
            ++front_of_camera_0;
        }

        // Point in the second camera coordinate system
        Vector3D P_cam1 = M_1 * P.homogeneous();
        if (P_cam1.z() > 0) {
            ++front_of_camera_1;
        }
    }

    std::cout << "Number of points in front of camera 0: " << front_of_camera_0 << "/" << points3D.size() << std::endl;
    std::cout << "Number of points in front of camera 1: " << front_of_camera_1 << "/" << points3D.size() << std::endl;
}

/*
 * This function is used for testing the four combinations of R and t.
 * The return integer is the number of reconstructed 3D points which
 * have positive z value in both cameras' frame (in front of the camera).
 */
int triangulation_test(
        const Matrix34& M_0, // input: K[I 0] matrix for the first camera.
        const Matrix34& M_1, // input: K[R, t] matrix, R and t is one of the four combinations.
        const std::vector<Vector2D> &points_0,  // input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  // input: 2D image points in the 2nd image.
        std::vector<Vector3D>& points3D // output: reconstructed 3d points from one combination.
) {
    int positive_z_count = 0;
    for (int i = 0; i < points_0.size(); i++) {
        Vector2D p0 = points_0[i];
        Vector2D p1 = points_1[i];
        Matrix44 A(4, 4, 0.0);
        Vector4D row0 = p0[0] * M_0.get_row(2) - M_0.get_row(0);
        Vector4D row1 = p0[1] * M_0.get_row(2) - M_0.get_row(1);
        Vector4D row2 = p1[0] * M_1.get_row(2) - M_1.get_row(0);
        Vector4D row3 = p1[1] * M_1.get_row(2) - M_1.get_row(1);
        A.set_row(0, row0);
        A.set_row(1, row1);
        A.set_row(2, row2);
        A.set_row(3, row3);

        // use SVD to get P
        Matrix44 U(4, 4, 0.0);
        Matrix44 D(4, 4, 0.0);
        Matrix44 V(4, 4, 0.0);
        svd_decompose(A, U, D, V);

        Vector4D P_homo = V.get_column(3);
        Vector3D P = P_homo.cartesian();
        points3D.push_back(P);

        /*
         * compute the coordinates of P under the second camera's frame.
         * note: this P2 is not the 3D point under the second camera's frame,
         * it is the homogeneous coordinates of the 2D point.
         */
        Vector3D P2 = M_1 * P.homogeneous();

        // check if the z values of P under both cameras are positive
        if (P[2] > 0 && P2[2] > 0) {
            positive_z_count += 1;
        }
    }
    return positive_z_count;
}

// This function performs triangulation using camera intrinsics (focal lengths, principal point, skew) and 2D image points, and outputs reconstructed 3D points, rotation matrix R, and translation vector t of the second camera.
bool Triangulation::triangulation(
        double fx, double fy,     // input: the focal lengths (same for both cameras)
        double cx, double cy,     // input: the principal point (same for both cameras)
        double s,                 // input: the skew factor (same for both cameras)
        const std::vector<Vector2D> &points_0,  // input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  // input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       // output: reconstructed 3D points
        Matrix33 &R,   // output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    // output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
// step 1: check inputs

// Ensure the number of 2D point pairs is consistent.
// Check for consistency in the number of points:
    if (points_0.size() != points_1.size()) {
        std::cerr << "Error: The number of points in the two images must be the same." << std::endl;
        return false;
    }

// Check for at least 8 points:
    if (points_0.size() < 8 || points_1.size() < 8) {
        std::cerr << "Error: At least 8 points are required for triangulation." << std::endl;
        return false;
    }

// Check for positive focal lengths:
    if (fx <= 0 || fy <= 0) {
        std::cerr << "Error: Focal lengths must be positive." << std::endl;
        return false;
    }

// Check the validity of the points:
    for (const Vector2D& p : points_0) {
        if (!std::isfinite(p[0]) || !std::isfinite(p[1])) {
            std::cerr << "Error: Invalid point in points_0." << std::endl;
            return false;
        }
    }

    for (const Vector2D& p : points_1) {
        if (!std::isfinite(p[0]) || !std::isfinite(p[1])) {
            std::cerr << "Error: Invalid point in points_1." << std::endl;
            return false;
        }
    }

// step 2: normalization of points

// Calculate the centroid of the input image points.
// Compute the mean distance of points from the centroid, then compute the scale factor so that the mean distance is sqrt(2). Normalize the points using the centroid and scale factor.

    // set up expected results: normalized points and T matrices
    std::vector<Vector2D> norm_points_0;
    std::vector<Vector2D> norm_points_1;
    Matrix33 T_0 = Matrix33::identity();
    Matrix33 T_1 = Matrix33::identity();

    // compute centroid of the input image points
    // This part of the code calculates the centroids (centroid_0 and centroid_1) for two sets of points (points_0 and points_1).
    Vector2D centroid_0(0, 0);
    Vector2D centroid_1(0, 0);
    for (const Vector2D& p: points_0) {
        centroid_0 += p;
    }
    centroid_0 /= points_0.size();

    for (const Vector2D& p: points_1) {
        centroid_1 += p;
    }
    centroid_1 /= points_1.size();

    // compute the scale factor for points_0 and points_1
    // This part of the code calculates the mean distance of each point from the centroid and then computes the scaling factor so that the average distance of the points is sqrt(2).
    double mean_dist_0 = 0;
    double mean_dist_1 = 0;

    for (const Vector2D& p: points_0) {
        mean_dist_0 += (p - centroid_0).length();
    }
    mean_dist_0 /= points_0.size();

    for (const Vector2D& p: points_1) {
        mean_dist_1 += (p - centroid_1).length();
    }
    mean_dist_1 /= points_1.size();

    double scale_0 = std::sqrt(2) / mean_dist_0;
    double scale_1 = std::sqrt(2) / mean_dist_1;

    // use centroids and scale factors to compute the two T matrices
    // Construct T matrices
    // This part of the code constructs two T matrices (T_0 and T_1) to translate and scale the points.
    T_0(0, 0) = scale_0;
    T_0(1, 1) = scale_0;
    T_0(0, 2) = - scale_0 * centroid_0.x();
    T_0(1, 2) = - scale_0 * centroid_0.y();

    T_1(0, 0) = scale_1;
    T_1(1, 1) = scale_1;
    T_1(0, 2) = - scale_1 * centroid_1.x();
    T_1(1, 2) = - scale_1 * centroid_1.y();

    // expand the scaled points to arrays
    // Normalize points
    // This part of the code multiplies the points by the T matrices to translate and scale them, resulting in normalized points.
    for (const Vector2D& p: points_0) {
        Vector2D normalized_p = T_0 * p.homogeneous();
        norm_points_0.push_back(normalized_p);
    }

    for (const Vector2D& p: points_1) {
        Vector2D normalized_p = T_1 * p.homogeneous();
        norm_points_1.push_back(normalized_p);
    }

// step 3: compute K matrix for both camera
// Use the input focal lengths, principal point, and skew to construct the intrinsic matrix K.
    Matrix33 K = Matrix33::identity();
    K(0, 0) = fx;
    K(1, 1) = fy;
    K(0, 1) = s;
    K(0, 2) = cx;
    K(1, 2) = cy;

// step 4: compute the fundamental matrix Fq using normalized points using SVD
// This part of the code describes how to calculate the fundamental matrix Fq using normalized points.
// The process includes constructing the equation matrix P, using singular value decomposition (SVD) to compute the fundamental matrix F, enforcing the rank-2 constraint, and then denormalizing it.

// Compute the fundamental matrix Fq using normalized points.
// Perform SVD to get Fq and enforce the rank-2 constraint.

    // compute the P matrix using normalized points
    // Construct the equation matrix P using normalized points.
    int m = norm_points_0.size();
    int n = 9;
    Matrix P(m, n, 0.0);
    for (int row = 0; row < norm_points_0.size(); row++) {
        Vector2D p = norm_points_0[row];
        Vector2D q = norm_points_1[row];
        double u = p[0];
        double v = p[1];
        double u1 = q[0];
        double v1 = q[1];

        P[row][0] = u * u1;
        P[row][1] = v * u1;
        P[row][2] = u1;
        P[row][3] = u * v1;
        P[row][4] = v * v1;
        P[row][5] = v1;
        P[row][6] = u;
        P[row][7] = v;
        P[row][8] = 1;
    }

    // compute U, D and V matrices and perform SVD to get scaled Fq: f
    // Compute Fq using SVD
    // The fundamental matrix F corresponds to the last column of V.
    Matrix U(m, m, 0.0);
    Matrix D(m, n, 0.0);
    Matrix V(n, n, 0.0);
    Matrix f(3, 3, 0.0);
    svd_decompose(P, U, D, V);
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            int index = row * 3 + col;
            f(row, col) = V.get_column(V.cols() - 1)[index];
        }
    }

    // perform SVD on f and do the constraint enforcement to get Fq
    // Perform SVD on f and set the last singular value to 0 to enforce the rank-2 constraint.
    Matrix33 U1(3, 3, 0.0);
    Matrix33 D1(3, 3, 0.0);
    Matrix33 V1(3, 3, 0.0);
    svd_decompose(f, U1, D1, V1);
    D1(2, 2) = 0.0; // Set the last singular value to 0 to enforce the rank-2 constraint
    Matrix33 Fq = U1 * D1 * inverse(V1);

// step 5: recover the expected fundamental matrix F by denormalization
// Denormalize Fq to get F:
    Matrix33 F = T_1.transpose() * Fq * T_0;

// step 6: compute the essential matrix E
// Calculate the essential matrix E using the formula E = K^T F K.
    Matrix33 E = K.transpose() * F * K;

// step 7: recover R and t // This part starts from step 2 of the assignment
// Perform SVD on the essential matrix E.
// Use the decomposition results to compute two possible rotation matrices and two possible translation vectors, for a total of four combinations.
// Reconstruct 3D points using each of the four combinations. Then, in step 8, calculate the number of 3D points with positive z values in both camera coordinate systems.

    // perform SVD on E to get U and V
    Matrix33 U2(3, 3, 0.0);
    Matrix33 D2(3, 3, 0.0);
    Matrix33 V2(3, 3, 0.0);
    svd_decompose(E, U2, D2, V2);

    // compute the W matrix
    // Construct the W matrix
    Matrix33 W(3, 3, 0.0);
    W.set_row(0, {0, -1, 0});
    W.set_row(1, {1, 0, 0});
    W.set_row(2, {0, 0, 1});

    // compute two potential R matrices
    // Calculate two possible rotation matrices
    Matrix33 SVD1 = U2 * W * inverse(V2);
    Matrix33 SVD2 = U2 * inverse(W) * inverse(V2);

    Matrix33 R1 = determinant(SVD1) * SVD1; // Ensure that the determinant of the rotation matrix R is 1 (orthogonal matrix).
    Matrix33 R2 = determinant(SVD2) * SVD2;

    // compute two potential t vectors
    // Calculate two possible translation vectors
    Vector3D t1 = U2.get_column(2);
    Vector3D t2 = - U2.get_column(2);

    // compute 4 combinations of R and t along with K
    // This step reconstructs 3D points using the four (R, t) combinations
    // and calculates the number of 3D points with positive z values in both camera coordinate systems.

    // Construct camera projection matrices
    Matrix34 M1(3, 4, 0.0);
    Matrix34 M2(3, 4, 0.0);
    Matrix34 M3(3, 4, 0.0);
    Matrix34 M4(3, 4, 0.0);

    M1.set_row(0, {R1(0, 0), R1(0, 1), R1(0, 2), t1[0]});
    M1.set_row(1, {R1(1, 0), R1(1, 1), R1(1, 2), t1[1]});
    M1.set_row(2, {R1(2, 0), R1(2, 1), R1(2, 2), t1[2]});
    M1 = K * M1;

    M2.set_row(0, {R1(0, 0), R1(0, 1), R1(0, 2), t2[0]});
    M2.set_row(1, {R1(1, 0), R1(1, 1), R1(1, 2), t2[1]});
    M2.set_row(2, {R1(2, 0), R1(2, 1), R1(2, 2), t2[2]});
    M2 = K * M2;

    M3.set_row(0, {R2(0, 0), R2(0, 1), R2(0, 2), t1[0]});
    M3.set_row(1, {R2(1, 0), R2(1, 1), R2(1, 2), t1[1]});
    M3.set_row(2, {R2(2, 0), R2(2, 1), R2(2, 2), t1[2]});
    M3 = K * M3;

    M4.set_row(0, {R2(0, 0), R2(0, 1), R2(0, 2), t2[0]});
    M4.set_row(1, {R2(1, 0), R2(1, 1), R2(1, 2), t2[1]});
    M4.set_row(2, {R2(2, 0), R2(2, 1), R2(2, 2), t2[2]});
    M4 = K * M4;

    // construct the M for the first camera
    Matrix34 Mf(3, 4, 0.0);
    Mf.set_row(0, {1, 0, 0, 0});
    Mf.set_row(1, {0, 1, 0, 0});
    Mf.set_row(2, {0, 0, 1, 0});
    Mf = K * Mf;

// step 8: reconstruct 3D points
// Reconstruct 3D points using each of the four combinations. Then, in step 8, calculate the number of 3D points with positive z values in both camera coordinate systems.
// Select the combination with the most points having positive z values as the final result, and return the reconstructed 3D points, rotation matrix R, and translation vector t.

    // Reconstruct 3D points and calculate the number of points with positive z values
    std::vector<Vector3D> points3D_1;
    std::vector<Vector3D> points3D_2;
    std::vector<Vector3D> points3D_3;
    std::vector<Vector3D> points3D_4;
    int positive_z_count_1 = triangulation_test(Mf, M1, points_0, points_1, points3D_1);
    int positive_z_count_2 = triangulation_test(Mf, M2, points_0, points_1, points3D_2);
    int positive_z_count_3 = triangulation_test(Mf, M3, points_0, points_1, points3D_3);
    int positive_z_count_4 = triangulation_test(Mf, M4, points_0, points_1, points3D_4);

    std::cout << "Positive z counts:" << std::endl;
    std::cout << "Combination 1: " << positive_z_count_1 << std::endl;
    std::cout << "Combination 2: " << positive_z_count_2 << std::endl;
    std::cout << "Combination 3: " << positive_z_count_3 << std::endl;
    std::cout << "Combination 4: " << positive_z_count_4 << std::endl;

    int final_score = std::max({positive_z_count_1, positive_z_count_2, positive_z_count_3, positive_z_count_4});

    std::vector<Vector3D> best_points3D;
    Matrix33 best_R;
    Vector3D best_t;
    Matrix34 best_M;

    // Select the best R and t
    if (final_score == positive_z_count_1) {
        best_points3D = points3D_1;
        best_R = R1;
        best_t = t1;
        best_M = M1;
        std::cout << "Using Combination 1 as the final result." << std::endl;
    } else if (final_score == positive_z_count_2) {
        best_points3D = points3D_2;
        best_R = R1;
        best_t = t2;
        best_M = M2;
        std::cout << "Using Combination 2 as the final result." << std::endl;
    } else if (final_score == positive_z_count_3) {
        best_points3D = points3D_3;
        best_R = R2;
        best_t = t1;
        best_M = M3;
        std::cout << "Using Combination 3 as the final result." << std::endl;
    } else {
        best_points3D = points3D_4;
        best_R = R2;
        best_t = t2;
        best_M = M4;
        std::cout << "Using Combination 4 as the final result." << std::endl;
    }

    points_3d = best_points3D;
    R = best_R;
    t = best_t;

    // step 9 evaluation
    Matrix34 M0 = K * Matrix34::identity();

    // Calculate the distribution of 3D points in both camera coordinate systems
    std::cout << " " << std::endl;
    std::cout << "Point distribution for Combination 1:" << std::endl;
    compute_point_distribution(points3D_1, M0, M1);

    std::cout << "Point distribution for Combination 2:" << std::endl;
    compute_point_distribution(points3D_2, M0, M2);

    std::cout << "Point distribution for Combination 3:" << std::endl;
    compute_point_distribution(points3D_3, M0, M3);

    std::cout << "Point distribution for Combination 4:" << std::endl;
    compute_point_distribution(points3D_4, M0, M4);

    // Calculate the distribution of 3D points in the best combination
    std::cout << " " << std::endl;
    std::cout << "Point distribution for the best combination:" << std::endl;
    compute_point_distribution(points_3d, M0, best_M);

    // Reprojection error calculation (best combination)
    double total_error = 0.0;
    for (int i = 0; i < points_3d.size(); ++i) {
        Vector4D P_homo(points_3d[i].x(), points_3d[i].y(), points_3d[i].z(), 1.0);

        Vector3D P_cam0 = M0 * P_homo; // Project the 3D point onto the first image plane
        Vector2D reprojected_p0(P_cam0.x() / P_cam0.z(), P_cam0.y() / P_cam0.z());

        Vector3D P_cam1 = best_M * P_homo; // Project the 3D point onto the second image plane
        Vector2D reprojected_p1(P_cam1.x() / P_cam1.z(), P_cam1.y() / P_cam1.z());

        double error0 = (reprojected_p0 - points_0[i]).length();
        double error1 = (reprojected_p1 - points_1[i]).length();

        total_error += error0 + error1;
    }
    double average_error = total_error / (2 * points_3d.size());
    std::cout << "Average reprojection error for the best combination: " << average_error << std::endl;

    return true;
}