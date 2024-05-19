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
void check_normalization(const std::vector<Vector2D>& points, const std::string& description);
Matrix construct_matrix_W(const std::vector<Vector2D>& points_0_normalized, const std::vector<Vector2D>& points_1_normalized);

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
    }


    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // Compute transformation matrices
    Matrix33 T0 = calculate_transformation_matrix(points_0);
    Matrix33 T1 = calculate_transformation_matrix(points_1);

    // Normalize points
    std::vector<Vector2D> points_0_normalized = normalize_points(points_0, T0);
    std::vector<Vector2D> points_1_normalized = normalize_points(points_1, T1);

    // Check normalization
    check_normalization(points_0_normalized, "Points 0");
    check_normalization(points_1_normalized, "Points 1");

    // construct matrix W
    Matrix W = construct_matrix_W(points_0_normalized, points_1_normalized);

    std::cout << "Matrix W has " << W.rows() << " rows and " << W.cols() << " columns" << std::endl;

    // Extract f matrix using SVD decomposition

    // First we need to extract the estimate, use the slides for this

    Matrix U(W.rows(), W.rows());
    Matrix D(W.rows(), W.cols());
    Matrix V(W.cols(), W.cols());


    std::cout << "Matrix U has " << U.rows() << " rows and " << U.cols() << " columns" << std::endl;
    std::cout << "Matrix D has " << D.rows() << " rows and " << D.cols() << " columns" << std::endl;
    std::cout << "Matrix V has " << V.rows() << " rows and " << V.cols() << " columns" << std::endl;

    svd_decompose(W, U, D, V);

    if (D.rows() >= 3) {
        D(2, 2) = 0.0;
    }

    Matrix F = U * D * V.transpose();

    std::cout << "Fundamental Matrix F:" << std::endl;
    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < F.cols(); ++j) {
            std::cout << F(i, j) << " ";
        }
        std::cout << std::endl;
    }




    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
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
        sumDistSquared += dx * dx + dy * dy;
    }
    double avgDist = std::sqrt(sumDistSquared / N);
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

void check_normalization(const std::vector<Vector2D>& points, const std::string& description) {
    double sumX = 0, sumY = 0;
    double sumDistSquared = 0;
    size_t N = points.size();

    for (const auto& p : points) {
        sumX += p.x();
        sumY += p.y();
        sumDistSquared += p.x() * p.x() + p.y() * p.y();
    }

    double centroidX = sumX / N;
    double centroidY = sumY / N;
    double avgDist = std::sqrt(sumDistSquared / N);

    std::cout << "Normalization check for " << description << ":" << std::endl;
    std::cout << " - Centroid: (" << centroidX << ", " << centroidY << ")" << std::endl;
    std::cout << " - Average Distance from Origin: " << avgDist << std::endl;

    // Check if the centroid is close enough to (0,0) and the average distance is close to sqrt(2)
    const double epsilon = 1e-6; // Tolerance for floating point comparison
    bool isCentroidCorrect = std::abs(centroidX) < epsilon && std::abs(centroidY) < epsilon;
    bool isScaleCorrect = std::abs(avgDist - std::sqrt(2)) < epsilon;

    std::cout << " - Centroid is " << (isCentroidCorrect ? "correct" : "incorrect") << std::endl;
    std::cout << " - Scale is " << (isScaleCorrect ? "correct" : "incorrect") << std::endl;
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