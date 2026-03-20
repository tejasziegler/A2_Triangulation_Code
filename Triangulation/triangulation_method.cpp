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

// ==================================== HELPER FUNCTIONS ============================================

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

Matrix34 construct_M(const Matrix33& K, const Matrix33& R, const Vector3D& t) {
    // Step 1: build [R | t] as a 3x4 matrix
    Matrix34 Rt(R(0,0), R(0,1), R(0,2), t.x(),
                R(1,0), R(1,1), R(1,2), t.y(),
                R(2,0), R(2,1), R(2,2), t.z());

    return K * Rt;
}

Matrix construct_A(const Vector2D& p0, const Vector2D& p1, const Matrix34& M0, const Matrix34& M1)
{
    Matrix A(4, 4, 0.0);

    double x  = p0.x();
    double y  = p0.y();
    double xp = p1.x();
    double yp = p1.y();

    // get_row returns a 4-element vector for each row of M
    A.set_row(0, x  * M0.get_row(2) - M0.get_row(0));
    A.set_row(1, y  * M0.get_row(2) - M0.get_row(1));
    A.set_row(2, xp * M1.get_row(2) - M1.get_row(0));
    A.set_row(3, yp * M1.get_row(2) - M1.get_row(1));

    return A;
}

Vector least_sq_from_A(const Matrix& A) {
    const int r = A.rows(); // number of coordinate pairs!
    const int c = 4;

    // Dimensions of matrices in SVD:
    Matrix U(r, r, 0.0);
    Matrix S(r, c, 0.0);
    Matrix V(c, c, 0.0);

    svd_decompose(A, U, S, V);

    // the last column (the solution vector P)
    // V.cols() - 1 is index 3
    Vector P = V.get_column(c-1);

    // From LiangLiang's helper notes! Much faster.
    return P;
}

std::vector<Vector3D> triangulate(const std::vector<Vector2D>& points_0,
                                   const std::vector<Vector2D>& points_1,
                                   const Matrix34& M0,
                                   const Matrix34& M1)
{
    std::vector<Vector3D> points;

    for (int i = 0; i < (int)points_0.size(); i++) {
        Matrix A = construct_A(points_0[i], points_1[i], M0, M1);
        Vector P_homo = least_sq_from_A(A);
        Vector3D point_3d(P_homo[0] / P_homo[3],
                          P_homo[1] / P_homo[3],
                          P_homo[2] / P_homo[3]);
        points.push_back(point_3d);
    }

    return points;
}

Vector2D reproject(const Vector3D& P3D, const Matrix34& M) {
    // Reprojection: M (3x4) * homogeneous point (4x1)
    Vector3D p_hat = (M * P3D.homogeneous());
    Vector2D p_cart = p_hat.cartesian();
    return p_cart;
}

double reproj_error(const Vector2D& p, const Vector3D& P3D, const Matrix34& M) {
    // Squared Euclidean distance in each image
    double sq_error = (p - reproject(P3D, M)).length2();
    return sqrt(sq_error); // returns RMSE!
}

// ==================================== LM OBJECTIVE ============================================

/// Levenberg-Marquardt method for non-linear least squares:

class ReprojectionObjective : public Objective_LM {
public:
    Matrix34 M0, M1;   /// projection matrices for camera 0 and camera 1
    Vector2D p0, p1;   /// observed 2D points in each image

    /// 4 residuals: (x,y) per camera
    /// 3 variables: X, Y, Z
    ReprojectionObjective() : Objective_LM(4, 3) {}

    int evaluate(const double *x, double *fvec) {
        Vector3D P3D(x[0], x[1], x[2]);

        Vector2D res0 = reproject(P3D, M0) - p0;
        Vector2D res1 = reproject(P3D, M1) - p1;

        fvec[0] = res0.x();
        fvec[1] = res0.y();
        fvec[2] = res1.x();
        fvec[3] = res1.y();

        return 0;
    }
};


// ================================= TRIANGULATION IMPLEMENTATION ======================================

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

    // --------------------- CHECK VALID INPUTS --------------------------------------------------------------------

    // 1. Enough point pairs (8-point algorithm requires at least 8)
    const int MIN_POINTS = 8;
    if (points_0.size() < MIN_POINTS || points_1.size() < MIN_POINTS) {
        std::cerr << "[Error] Not enough point pairs: got " << points_0.size()
                  << " (image 0) and " << points_1.size()
                  << " (image 1), but at least " << MIN_POINTS << " are required." << std::flush;
        return false;
    }

    // 2. Equal number of points in both images
    if (points_0.size() != points_1.size()) {
        std::cerr << "[Error] Point count mismatch: image 0 has " << points_0.size()
                  << " points, image 1 has " << points_1.size() << " points." << std::flush;
        return false;
    }

    // All 2D point coordinates must be finite real numbers
    //        Vector2D inherits operator[] from Vector, giving access to x (index 0) and y (index 1).
    for (std::size_t i = 0; i < points_0.size(); ++i) {
        if (!std::isfinite(points_0[i][0]) || !std::isfinite(points_0[i][1])) {
            std::cerr << "[Error] Non-finite coordinate in image 0 at index " << i
                      << ": (" << points_0[i][0] << ", " << points_0[i][1] << ")" << std::flush;
            return false;
        }
        if (!std::isfinite(points_1[i][0]) || !std::isfinite(points_1[i][1])) {
            std::cerr << "[Error] Non-finite coordinate in image 1 at index " << i
                      << ": (" << points_1[i][0] << ", " << points_1[i][1] << ")" << std::flush;
            return false;
        }
    }

    // Camera intrinsic parameters must be valid -----------------------
    //        Focal lengths must be strictly positive finite values.
    //        Principal point and skew must be finite (any real value is geometrically valid).
    if (!std::isfinite(fx) || fx <= 0.0) {
        std::cerr << "[Error] Invalid focal length fx = " << fx
                  << ". Must be a positive finite number." << std::flush;
        return false;
    }
    if (!std::isfinite(fy) || fy <= 0.0) {
        std::cerr << "[Error] Invalid focal length fy = " << fy
                  << ". Must be a positive finite number." << std::flush;
        return false;
    }
    if (!std::isfinite(cx)) {
        std::cerr << "[Error] Invalid principal point cx = " << cx
                  << ". Must be a finite number." << std::flush;
        return false;
    }
    if (!std::isfinite(cy)) {
        std::cerr << "[Error] Invalid principal point cy = " << cy
                  << ". Must be a finite number." << std::flush;
        return false;
    }
    if (!std::isfinite(s)) {
        std::cerr << "[Error] Invalid skew factor s = " << s
                  << ". Must be a finite number." << std::flush;
        return false;
    }

    std::cout << "\n[1/6] Input validation passed: n=" << points_0.size() << " point pairs.\n";

    // --------------------- NORMALIZATION ----------------------------------------------------------------

    // Getting T, T', q, and q' using previously defined helper function
    std::vector<Vector3D> q0, q1;
    Matrix33 T  = normalize_points(points_0, q0);
    Matrix33 Tprime = normalize_points(points_1, q1);

    std::cout << "[2/6] Points normalized.\n";

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
    std::cout << "[3/6] Computed fundamental matrix F" << std::endl;

    // --------------------- FIND R AND t FROM F -------------------------------------------------------


    // From F to E
    // 1. Construct the Intrinsic Matrix K'
    Matrix33 K_1(fx,   s,  cx,
                 0,  fy,  cy,
                 0,   0,   1 );
    // 2. Construct the Intrinsic Matrix K
    Matrix33 K_0(fx,   s,  cx,
                 0,  fy,  cy,
                 0,   0,   1 );
    // 3. Derived and decompose the E

    Matrix33 E = K_0.transpose() * F * K_0;

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

    // Translation options: third column of U
    Vector3D tA = Vector3D(U1.get_column(2)); // get_column is 0-indexed, so 2 is the 3rd column
    Vector3D tB = -tA;

    std::cout << "[4/6] Essential matrix E and R/t candidates extracted.\n";


    // --------------------- TRIANGULATION -------------------------------------------------------------

    // camera 0 — fixed reference
    Matrix33 R0 = Matrix::identity(3, 3, 1.0);
    Vector3D t0(0, 0, 0);
    Matrix34 M0 = construct_M(K_0, R0, t0);

    // 4 candidate combinations
    std::vector<Matrix33> R_candidates = {RA, RB, RA, RB};
    std::vector<Vector3D> t_candidates = {tA, tA, tB, tB};

    // store the 4 resulting point sets
    std::vector<std::vector<Vector3D>> all_points_3d(4);

    for (int i = 0; i < 4; i++) {
        Matrix34 M1 = construct_M(K_1, R_candidates[i], t_candidates[i]);
        all_points_3d[i] = triangulate(points_0, points_1, M0, M1);
        std::cout << "[5/6] Candidate " << i << ": triangulated "
              << all_points_3d[i].size() << " points.\n";
    }

    // pick the set with the most points with positive z
    int best = 0;
    int best_count = 0;
    for (int i = 0; i < 4; i++) {
        int count = 0;
        for (const auto& p : all_points_3d[i]) {
            // depth in camera 0: just z (since R0=I, t0=0 in this case)
            Vector3D p_cam0 = R0 * p + t0;
            bool in_front_of_cam0 = p_cam0.z() > 0;

            // depth in camera 1: transform point into camera 1 space first
            Vector3D p_cam1 = R_candidates[i] * p + t_candidates[i];
            bool in_front_of_cam1 = p_cam1.z() > 0;

            if (in_front_of_cam0 && in_front_of_cam1) count++;
        }
        std::cout << "      Candidate " << i << ": " << count
              << " / " << all_points_3d[i].size() << " points in front of both cameras.\n";
        if (count > best_count) { best_count = count; best = i; }
    }

    // best solution
    R = R_candidates[best];
    t = t_candidates[best];
    points_3d = all_points_3d[best];

    std::cout << "[6/6] Best: candidate " << best
          << " with " << best_count << " valid points.\n";

    // Compute reprojection error
    Matrix34 M0_final = construct_M(K_0, R0, t0);
    Matrix34 M1_final = construct_M(K_1, R, t);

    double total_error = 0.0;
    for (size_t i = 0; i < points_3d.size(); ++i) {
        total_error += reproj_error(points_0[i], points_3d[i], M0_final) + reproj_error(points_1[i], points_3d[i], M1_final);
    }
    double rmse = total_error / (2*points_3d.size());
    // division by total number of 2D points
    std::cout << "     Root Mean Squared Reprojection Error: " << rmse << " squared pixels\n";

    // --------------------- LM IMPROVEMENT --------------------------------------------------

    int refined = 0;

    for (int i = 0; i < (int)points_3d.size(); i++) {

        /// initialize the objective function
        /// 1st argument is the number of functions, 2nd the number of variables
        ReprojectionObjective obj;
        obj.M0 = M0_final;
        obj.M1 = M1_final;
        obj.p0 = points_0[i];
        obj.p1 = points_1[i];

        /// create an instance of the Levenberg-Marquardt (LM for short) optimizer
        Optimizer_LM lm;

        Optimizer_LM::Parameters params;
        params.ftol    = 1e-10;
        params.xtol    = 1e-10;
        params.gtol    = 1e-10;
        params.maxcall = 1000;
        params.epsilon = 1e-8;

        /// initialize the variables from the linear triangulation result
        std::vector<double> x = { points_3d[i].x(),
                                   points_3d[i].y(),
                                   points_3d[i].z() };

        /// optimize (i.e., minimizing the reprojection error)
        bool status = lm.optimize(&obj, x, &params);

        /// retrieve the refined result
        if (status) {
            points_3d[i] = Vector3D(x[0], x[1], x[2]);
            refined++;
        }
    }

    std::cout << "      LM refinement: " << refined << " / "
              << points_3d.size() << " points refined.\n";

    // Recompute MSE after refinement
    double total_error_refined = 0.0;
    for (size_t i = 0; i < points_3d.size(); ++i) {
        total_error_refined += reproj_error(points_0[i], points_3d[i], M0_final)
                             + reproj_error(points_1[i], points_3d[i], M1_final);
    }
    std::cout << "      RMSE after LM: "
              << total_error_refined / (2 * points_3d.size()) << " squared pixels\n";
    // division by total number of 2D points

    return points_3d.size() > 0;
    }

