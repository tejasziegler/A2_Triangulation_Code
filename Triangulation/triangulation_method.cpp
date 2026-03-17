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

Matrix33 normalize_points(const std::vector<Vector2D>& pts, std::vector<Vector3D>& normalized) {
    int n = pts.size();

    double mean_x = 0, mean_y = 0;
    for (auto& p : pts) { mean_x += p[0]; mean_y += p[1]; }
    mean_x /= n;
    mean_y /= n;

    double mean_dist = 0;
    for (auto& p : pts)
        mean_dist += std::sqrt(std::pow(p[0] - mean_x, 2) + std::pow(p[1] - mean_y, 2));
    mean_dist /= n;

    double scale = std::sqrt(2.0) / mean_dist;
    Matrix33 T(scale, 0,     -scale * mean_x,
               0,     scale, -scale * mean_y,
               0,     0,      1);

    normalized.resize(n);
    for (int i = 0; i < n; ++i)
        normalized[i] = T * pts[i].homogeneous();

    return T;
}

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
) const {
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       'triangulation()'.

    std::cout << "\nTODO: implement the 'triangulation()' function in the file 'Triangulation/triangulation_method.cpp'\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tSimilar to the first assignment, basic linear algebra data structures and functions are provided in\n"
                 "\tthe following files:\n"
                 "\t    - Triangulation/matrix.h: handles matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h: manages vectors of arbitrary sizes and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h: contains functions for determinant, inverse, SVD, linear least-squares...\n"
                 "\tFor more details about these data structures and a complete list of related functions, please\n"
                 "\trefer to the header files mentioned above.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tFor your final submission, adhere to the following guidelines:\n"
                 "\t    - submit ONLY the 'Triangulation/triangulation_method.cpp' file.\n"
                 "\t    - remove ALL unrelated test code, debugging code, and comments.\n"
                 "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" << std::flush;

    // /// Below are a few examples showing some useful data structures and APIs.
    //
    // /// define a 2D vector/point
    // Vector2D b(1.1, 2.2);
    //
    // /// define a 3D vector/point
    // Vector3D a(1.1, 2.2, 3.3);
    //
    // /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    // Vector2D p = a.cartesian();
    //
    // /// get the Homogeneous coordinates of p
    // Vector3D q = p.homogeneous();
    //
    // /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    // Matrix33 A;
    //
    // /// define and initialize a 3 by 3 matrix
    // Matrix33 T(1.1, 2.2, 3.3,
    //            0, 2.2, 3.3,
    //            0, 0, 1);
    //
    // /// define and initialize a 3 by 4 matrix
    // Matrix34 M(1.1, 2.2, 3.3, 0,
    //            0, 2.2, 3.3, 1,
    //            0, 0, 1, 1);
    //
    // /// set first row by a vector
    // M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));
    //
    // /// set second column by a vector
    // M.set_column(1, Vector3D(5.5, 5.5, 5.5));
    //
    // /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    // Matrix W(15, 9, 0.0);
    // /// set the first row by a 9-dimensional vector
    // W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>
    //
    // /// get the number of rows.
    // int num_rows = W.rows();
    //
    // /// get the number of columns.
    // int num_cols = W.cols();
    //
    // /// get the the element at row 1 and column 2
    // double value = W(1, 2);
    //
    // /// get the last column of a matrix
    // Vector last_column = W.get_column(W.cols() - 1);
    //
    // /// define a 3 by 3 identity matrix
    // Matrix33 I = Matrix::identity(3, 3, 1.0);
    //
    // /// matrix-vector product
    // Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4
    //
    // ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).

    // 1. Enough point pairs (8-point algorithm requires at least 8)
    const int MIN_POINTS = 8;
    if (points_0.size() < MIN_POINTS || points_1.size() < MIN_POINTS) {
        std::cerr << "[Error] Not enough point pairs: got " << points_0.size()
                  << " (image 0) and " << points_1.size()
                  << " (image 1), but at least " << MIN_POINTS << " are required." << std::endl;
        return false;
    }

    // 2. Equal number of points in both images
    if (points_0.size() != points_1.size()) {
        std::cerr << "[Error] Point count mismatch: image 0 has " << points_0.size()
                  << " points, image 1 has " << points_1.size() << " points." << std::endl;
        return false;
    }

    // All 2D point coordinates must be finite real numbers
    //        Vector2D inherits operator[] from Vector, giving access to x (index 0) and y (index 1).
    for (std::size_t i = 0; i < points_0.size(); ++i) {
        if (!std::isfinite(points_0[i][0]) || !std::isfinite(points_0[i][1])) {
            std::cerr << "[Error] Non-finite coordinate in image 0 at index " << i
                      << ": (" << points_0[i][0] << ", " << points_0[i][1] << ")" << std::endl;
            return false;
        }
        if (!std::isfinite(points_1[i][0]) || !std::isfinite(points_1[i][1])) {
            std::cerr << "[Error] Non-finite coordinate in image 1 at index " << i
                      << ": (" << points_1[i][0] << ", " << points_1[i][1] << ")" << std::endl;
            return false;
        }
    }

    // Camera intrinsic parameters must be valid -----------------------
    //        Focal lengths must be strictly positive finite values.
    //        Principal point and skew must be finite (any real value is geometrically valid).
    if (!std::isfinite(fx) || fx <= 0.0) {
        std::cerr << "[Error] Invalid focal length fx = " << fx
                  << ". Must be a positive finite number." << std::endl;
        return false;
    }
    if (!std::isfinite(fy) || fy <= 0.0) {
        std::cerr << "[Error] Invalid focal length fy = " << fy
                  << ". Must be a positive finite number." << std::endl;
        return false;
    }
    if (!std::isfinite(cx)) {
        std::cerr << "[Error] Invalid principal point cx = " << cx
                  << ". Must be a finite number." << std::endl;
        return false;
    }
    if (!std::isfinite(cy)) {
        std::cerr << "[Error] Invalid principal point cy = " << cy
                  << ". Must be a finite number." << std::endl;
        return false;
    }
    if (!std::isfinite(s)) {
        std::cerr << "[Error] Invalid skew factor s = " << s
                  << ". Must be a finite number." << std::endl;
        return false;
    }




    // Getting T, T', q, and q' using previously defined helper function
    std::vector<Vector3D> q0, q1;
    Matrix33 T  = normalize_points(points_0, q0);
    Matrix33 Tprime = normalize_points(points_1, q1);

    // Build W
    int n = q0.size();
    Matrix W(n, 9, 0.0);
    for (int i = 0; i < n; ++i) {
        double u  = q0[i][0],  v  = q0[i][1];
        double u1 = q1[i][0],  v1 = q1[i][1];
        W.set_row(i, {u*u1, v*u1, u1, u*v1, v*v1, v1, u, v, 1.0});
    }

    // Solve Wf = 0 with SVD
    Matrix U(n, n, 0.0), S(n, 9, 0.0), V(9, 9, 0.0);
    svd_decompose(W, U, S, V);
    Vector f = V.get_column(V.cols() - 1);
    Matrix33 Fq(f[0], f[1], f[2],
                f[3], f[4], f[5],
                f[6], f[7], f[8]);

    // Enforce rank-2
    Matrix33 U2, V2;
    Matrix S2(3, 3, 0.0);
    svd_decompose(Fq, U2, S2, V2);
    S2(2, 2) = 0.0;
    Fq = U2 * S2 * transpose(V2);

    // De-normalize
    Matrix33 F = transpose(Tprime) * Fq * T;
    std::cout << "F (final fundamental matrix):\n" << F << std::endl;

    // Sanity check: for a correct F, p'^T * F * p should be near 0 for all point pairs
    // Sanity check
    std::cout << "Epipolar constraint check (should be near 0 for all pairs):" << std::endl;
    for (int i = 0; i < std::min(n, 5); ++i) {
        Vector3D p  = points_0[i].homogeneous();
        Vector3D p1 = points_1[i].homogeneous();
        double constraint = dot(p1, F * p);
        std::cout << "  pair " << i << ": " << constraint << std::endl;
    }  // <-- for loop ends HERE

    // 1. Extract R an t from F
    
    // From F to E
    // 1. Construct the Intrinsic Matrix K'
    Matrix33 K_1(fx,   s,  cx,
                 0,  fy,  cy,
                 0,   0,   1 );
    // 2. Construct the Intrinsic Matrix K
    Matrix33 K_0(1,   0,   0,
                 0,   1,   0,
                 0,   0,   1 );
    // 3. Derived and decompose the E

    Matrix33 E = K_1.transpose() * F * K_1;

    int r = E.rows();
    int c = E.cols();

    // Dimensions of matrices in SVD:
    Matrix U1(r, r, 0.0);
    Matrix S1(r, c, 0.0);
    Matrix V1(c, c, 0.0);

    svd_decompose(E, U1, S1, V1);

    // 4. Extract R and t candidates

    // Step A: Define the Helper Matrix W
    Matrix33 W1(0, -1, 0,
               1,  0, 0,
               0,  0, 1);

    // Step B: Extract the 4 Candidates
    Matrix33 RA = Matrix33(U1 * W1 * V1.transpose());
    Matrix33 RB = Matrix33(U1 * W1.transpose() * V1.transpose());

    // Ensure RA and RB are valid rotation matrices (determinant must be +1)
    if (determinant(RA) < 0) {
        RA = RA * -1.0;
    }
    if (determinant(RB) < 0) {
        RB = RB * -1.0;
    }

    std::cout << "RA " << RA << std::endl;
    std::cout << "RB " << RB << std::endl;  

    // Translation options: third column of U
    Vector3D tA = Vector3D(U.get_column(2)); // get_column is 0-indexed, so 2 is the 3rd column
    Vector3D tB = -tA;
    std:: cout << "t " << tA << std::endl;
    std:: cout << "t " << tB << std::endl;
        // TODO: Estimate relative pose of two views. This can be subdivided into
        //      - estimate the fundamental matrix F;
        //      - compute the essential matrix E;
        //      - recover rotation R and t.

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



