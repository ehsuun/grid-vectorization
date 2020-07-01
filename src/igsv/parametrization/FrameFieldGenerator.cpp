#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <igl/LinSpaced.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/per_vertex_attribute_smoothing.h>
#include <igl/slice_into.h>
#include <igl/speye.h>
#include <igsv/parametrization/FrameFieldGenerator.h>
#include <iostream>

namespace IGSV {

  using complexd     = std::complex<double>;
  using SparseMatrix = Eigen::SparseMatrix<double>;
  using Tripletd     = Eigen::Triplet<double>;
  using VecXcd       = Eigen::VectorXcd;
  using VecXd        = Eigen::VectorXd;
  using VecXi        = Eigen::VectorXi;

  // ===========================================================================

  void FrameFieldGenerator::init() {
    V_constraints.resize(_V.rows());
    int c = 0;
    for (auto v : _C_nearestVertices)
      V_constraints[v].push_back(c++);
  }

  // ===========================================================================

  bool FrameFieldGenerator::compute() {
    VecXd X_;
    if (!solve(X_))
      return false;
    return extract(X_);
  }

  // ===========================================================================

  void FrameFieldGenerator::extract_roots(const complexd& c0, //
                                          const complexd& c2, //
                                          double& X0_arg,     //
                                          double& X0_mag,     //
                                          double& X1_arg,     //
                                          double& X1_mag) {
    /*
      conversion formulas from [Bessmeltsev & Solomon 2018]
      "Vectorization of Line Drawings via PolyVector Fields"
      https://arxiv.org/abs/1801.01922
  */
    const complexd a(+4.0, 0.);
    const complexd b(-0.5, 0.);
    const complexd u2 = b * (c2 + std::sqrt(c2 * c2 - a * c0));
    const complexd v2 = b * (c2 - std::sqrt(c2 * c2 - a * c0));

    X0_arg = 0.5 * std::arg(u2);
    X0_mag = std::sqrt(std::abs(u2));
    X1_arg = 0.5 * std::arg(v2);
    X1_mag = std::sqrt(std::abs(v2));

    // make sure we take the correct pair
    if (std::abs(X0_arg - X1_arg) > 0.5 * M_PI) {
      X1_arg += (X1_arg > 0) ? -M_PI : M_PI;
    }
  }

  // ===========================================================================

  bool FrameFieldGenerator::build_energy_smooth(SparseMatrix& A, VecXd& B) {

    const int nv = _V.rows();
    const int nf = _F.rows();

    //// triangle areas
    VecXd W(nf);
    W.fill(1.0);
    igl::doublearea(_V, _F, W);

    //// gradient
    SparseMatrix G(nf * 3, nv);
    igl::grad(_V, _F, G);

    //// weighted Laplacian
    SparseMatrix L = 0.5 * G.transpose() * W.replicate<3, 1>().asDiagonal() * G;

    //// energy matrix
    A.resize(4 * nv, 4 * nv);
    for (int k = 0; k < 4; k++) {
      const VecXi idx = igl::LinSpaced<VecXi>(nv, k * nv, (k + 1) * nv - 1);
      igl::slice_into(L, idx, idx, A);
    }

    //// right-hand side
    B.setZero(4 * nv, 1);

    return true;
  }

  // ===========================================================================

  bool FrameFieldGenerator::build_energy_align(const VecXcd& tau,    //
                                               const VecXd& weights, //
                                               const VecXi& nearest, //
                                               int nv,               //
                                               SparseMatrix& A,      //
                                               VecXd& B,             //
                                               bool grad_equal_zero) {

    const int nc = tau.rows();

    assert(weights.rows() == nc);
    assert(nearest.rows() == nc);

    const auto x = tau.real();
    const auto y = tau.imag();

    const auto xy0 = 2 * x.array() * y.array();
    const auto xy1 = x.array().pow(2) - y.array().pow(2);
    const auto xy2 = x.array().pow(4) + 2 * x.array().pow(2) * y.array().pow(2) + y.array().pow(2);

    std::vector<Tripletd> ijv;
    ijv.reserve(16 * nc);

    for (int k = 0; k < nc; ++k) {
      const int v = nearest(k);
      if (v < 0 || v >= nv)
        return false;

      // rows
      const int r0 = k + 0 * nc;
      const int r1 = k + 1 * nc;

      // cols
      const int x0 = v + 0 * nv;
      const int y0 = v + 1 * nv;
      const int x2 = v + 2 * nv;
      const int y2 = v + 3 * nv;

      // if (std::isnan(xy0(k))) printf("WARNING : xy0(k) = NaN (k=%d)\n",k);
      // if (std::isnan(xy1(k))) printf("WARNING : xy1(k) = NaN (k=%d)\n",k);
      // if (std::isnan(xy2(k))) printf("WARNING : xy2(k) = NaN (k=%d)\n",k);

      ijv.emplace_back(r0, x0, +1.0);
      ijv.emplace_back(r0, y0, 0.0);
      ijv.emplace_back(r0, x2, +xy1(k));
      ijv.emplace_back(r0, y2, -xy0(k));

      ijv.emplace_back(r1, x0, 0.0);
      ijv.emplace_back(r1, y0, +1.0);
      ijv.emplace_back(r1, x2, +xy0(k));
      ijv.emplace_back(r1, y2, +xy1(k));

      if (grad_equal_zero) {
        const int r2 = k + 2 * nc;
        const int r3 = k + 3 * nc;
        ijv.emplace_back(r2, x0, +xy1(k));
        ijv.emplace_back(r2, y0, +xy0(k));
        ijv.emplace_back(r2, x2, +xy2(k));
        ijv.emplace_back(r2, y2, 0.0);

        ijv.emplace_back(r3, x0, -xy0(k));
        ijv.emplace_back(r3, y0, +xy1(k));
        ijv.emplace_back(r3, x2, 0.0);
        ijv.emplace_back(r3, y2, +xy2(k));
      }
    }

    const int nEquationsPerConstraint = grad_equal_zero ? 4 : 2;

    // Lhs
    SparseMatrix Lhs;
    Lhs.reserve(16 * nc);
    Lhs.resize(nEquationsPerConstraint * nc, 4 * nv);
    Lhs.setFromTriplets(ijv.begin(), ijv.end());

    // Rhs
    VecXd Rhs;
    Rhs.resize(nEquationsPerConstraint * nc, 1);
    Rhs.segment(0 * nc, nc) = x.array().pow(4) - 6 * y.array().pow(2) * x.array().pow(2) + y.array().pow(4);
    Rhs.segment(1 * nc, nc) = 4 * x.array().pow(3) * y.array() - 4 * x.array() * y.array().pow(3);

    if (grad_equal_zero) {
      Rhs.segment(2 * nc, nc) = x.array().pow(6) + x.array().pow(4) * y.array().pow(2) -
                                x.array().pow(2) * y.array().pow(4) - y.array().pow(6);
      Rhs.segment(3 * nc, nc) =
          2 * x.array().pow(5) * y.array() + 4 * x.array().pow(3) * y.array().pow(3) + 2 * x.array() * y.array().pow(5);
    }
    Rhs.array() *= -1.; // flip

    // Energy terms with weights
    if (grad_equal_zero) {
      A = Lhs.transpose() * weights.replicate<4, 1>().asDiagonal() * Lhs;
      B = Lhs.transpose() * weights.replicate<4, 1>().asDiagonal() * Rhs;
    } else {
      A = Lhs.transpose() * weights.replicate<2, 1>().asDiagonal() * Lhs;
      B = Lhs.transpose() * weights.replicate<2, 1>().asDiagonal() * Rhs;
    }

    return true;
  }

  // ===========================================================================

  bool FrameFieldGenerator::solve(SparseMatrix& A_smooth,  //
                                  SparseMatrix& A_tangent, //
                                  SparseMatrix& A_ortho,   //
                                  VecXd& B_smooth,         //
                                  VecXd& B_tangent,        //
                                  VecXd& B_ortho,          //
                                  VecXd& X_) {

    const int nv = _V.rows();

    VecXd X_prev;

    if (!build_energy_smooth(A_smooth, B_smooth))
      return false;

    if (!build_energy_align(_C_tangents,        //
                            _C_tangentWeights,  //
                            _C_nearestVertices, //
                            nv,                 //
                            A_tangent,          //
                            B_tangent,          //
                            true))
      return false;

    const VecXcd C_normals = _C_tangents * complexd(0.0, 1.0); // normal = R90 * tangent

    // init normal weights
    Eigen::VectorXd C_normalWeights;
    C_normalWeights.resize(_C_tangentWeights.rows(), 1);
    if (_use_weighted_reg)
      C_normalWeights.fill(0.);
    else
      C_normalWeights.fill(1.);

#define MAX_ITERS 100
#define MAG_RATIO_LIMIT 0.25
#define N_DEGENERATE_RATIO_LIMIT 0.005

    const double ANGLE_DIF_LIMIT = this->_minangle / 180. * M_PI;

    double uangle, umag; // first frame vector
    double vangle, vmag; // second frame vector
    double angle_dif;    // angle difference
    double mag_ratio;    // magnitude ratio

    bool are_there_degenerate_frames = true;

    int n_degenerate_prev;

    _n_iters = 0;
    while (are_there_degenerate_frames && _n_iters++ < MAX_ITERS) {

      // normal alignment energy with the current weights
      if (!build_energy_align(C_normals,          //
                              C_normalWeights,    //
                              _C_nearestVertices, //
                              nv,                 //
                              A_ortho,            //
                              B_ortho,            //
                              true))
        return false;

      // small regularization
      SparseMatrix A_reg;
      igl::speye(4 * nv, 4 * nv, A_reg);

      SparseMatrix A_ = _wsmooth * A_smooth + _wtangent * A_tangent + _wortho * A_ortho + 0.0001 * A_reg;
      VecXd B_        = _wsmooth * B_smooth + _wtangent * B_tangent + _wortho * B_ortho;

      // solve via Cholesky factorization
      Eigen::SimplicialLDLT<SparseMatrix> chol_solver(A_);

      if (chol_solver.info() == Eigen::Success) {

        X_ = chol_solver.solve(B_);
        if (chol_solver.info() != Eigen::Success) {
          return false;
        }

      } else {
        // fallback to QR solver
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr_solver;
        A_.makeCompressed();
        qr_solver.compute(A_);

        if (qr_solver.info() == Eigen::Success) {
          X_ = qr_solver.solve(B_);

          if (qr_solver.info() != Eigen::Success) {
            return false;
          }
        } else {
          return false;
        }
      }

      // check for degenerate triangles!
      are_there_degenerate_frames = false;
      _n_degenerate               = 0;

      for (int v = 0; v < nv; ++v) {

        extract_roots(complexd(X_[v + 0 * nv], X_[v + 1 * nv]), ///
                      complexd(X_[v + 2 * nv], X_[v + 3 * nv]), //
                      uangle, umag, vangle, vmag);

        angle_dif = std::abs(uangle - vangle);
        while (angle_dif > +M_PI)
          angle_dif -= 2 * M_PI;
        while (angle_dif < -M_PI)
          angle_dif += 2 * M_PI;

        if (umag < 1e-4 && vmag < 1e-4)
          mag_ratio = 0.0;
        else
          mag_ratio = (umag > vmag) ? vmag / umag : umag / vmag;

        if (angle_dif < ANGLE_DIF_LIMIT || mag_ratio < MAG_RATIO_LIMIT) {
          for (auto c : V_constraints[v])
            C_normalWeights(c) += 1.0;
          are_there_degenerate_frames = true;
          _n_degenerate++;
        }
      }

      if (_n_iters > 1) {
        // check how many degenerate frames have been removed
        const double n_removed_ratio    = (double)(n_degenerate_prev - _n_degenerate) / (double)(n_degenerate_prev);
        const double n_degenerate_ratio = (double)_n_degenerate / (double)nv;
        // std::cout << "    " << n_degenerate << " / " << nv << " degenerate frames    (removed "
        //           << std::round<int>(n_removed_ratio * 100.) << "%)\n\n";

        if (n_removed_ratio < 0.0) {
          // the solution actually got worse, revert and stop
          X_ = X_prev;
          break;
        }

        if (n_degenerate_ratio < N_DEGENERATE_RATIO_LIMIT)
          break;

      } else {
        // std::cout << "    " << n_degenerate << " / " << nv << " degenerate frames" << std::endl << std::endl;
      }

      X_prev = X_;

      n_degenerate_prev = _n_degenerate;

      if (!_use_weighted_reg)
        break;

    } // end (while there are degenerate frames)

    return true;
  }

  // ===========================================================================

  bool FrameFieldGenerator::solve(VecXd& X_) {

    SparseMatrix A_smooth, A_tangent, A_ortho;
    VecXd B_smooth, B_tangent, B_ortho;

    return solve(A_smooth, A_tangent, A_ortho, B_smooth, B_tangent, B_ortho, X_);
  }

  // ===========================================================================

  bool FrameFieldGenerator::extract(const VecXd& X) {

    const int nv = _V.rows();
    const int nf = _F.rows();

    // store per-vertex frame field
    _X0_arg_vtx.resize(nv, 1);
    _X1_arg_vtx.resize(nv, 1);
    _X0_mag_vtx.resize(nv, 1);
    _X1_mag_vtx.resize(nv, 1);

    for (int v = 0; v < nv; ++v) {

      const complexd c0(X[v + 0 * nv], X[v + 1 * nv]);
      const complexd c2(X[v + 2 * nv], X[v + 3 * nv]);

      double x0_arg, x0_mag, x1_arg, x1_mag;
      extract_roots(c0, c2, x0_arg, x0_mag, x1_arg, x1_mag);

      _X0_arg_vtx(v) = x0_arg;
      _X0_mag_vtx(v) = x0_mag;
      _X1_arg_vtx(v) = x1_arg;
      _X1_mag_vtx(v) = x1_mag;
    }

    // store per-face frame field
    _X0_arg.resize(nf, 1);
    _X1_arg.resize(nf, 1);
    _X0_mag.resize(nf, 1);
    _X1_mag.resize(nf, 1);

    for (int f = 0; f < nf; ++f) {
      double c0_real = 0;
      double c0_imag = 0;
      double c2_real = 0;
      double c2_imag = 0;

      for (int k = 0; k < 3; ++k) {
        const int i = _F(f, k);
        c0_real += X[i + 0 * nv];
        c0_imag += X[i + 1 * nv];
        c2_real += X[i + 2 * nv];
        c2_imag += X[i + 3 * nv];
      }

      c0_real /= 3.0;
      c0_imag /= 3.0;
      c2_real /= 3.0;
      c2_imag /= 3.0;

      // extract roots
      const complexd c0(c0_real, c0_imag);
      const complexd c2(c2_real, c2_imag);

      double x0_arg, x0_mag, x1_arg, x1_mag;
      extract_roots(c0, c2, x0_arg, x0_mag, x1_arg, x1_mag);

      _X0_arg(f) = x0_arg;
      _X0_mag(f) = x0_mag;
      _X1_arg(f) = x1_arg;
      _X1_mag(f) = x1_mag;
    }
    return true;
  }

  // ===========================================================================

} // namespace IGSV
