//=============================================================================
//
//  CLASS : PoissonElement
//          SnappingElement
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <cfloat>
#include <cmath>
#include <complex>
#include <iostream>

//== NAMESPACES ===============================================================

namespace IGSV {

  //===========================================================================
  //== CLASS DEFINITION =======================================================
  //===========================================================================
  // Poisson Energy Element
  // element for 1/2 * ||Df(F) - I||_alpha_^2_{Fro} = 1/(4A) * F^t * G * U - I
  // with
  //  U = (u_0, v_0; u_1, v_1; u_2, v_2) : 3x2 matrix of unknown uv coordinates
  //  G = Rot90*(e_0, e_1, e_2)          : 2x3 gradient matrix
  //  F = (x_u, x_v; y_u, y_v)           : 2x2 frame matrix
  //  I = (1, 0; 0, 1)                   : 2x2 identity matrix
  class PoissonElement {

  public:
    // define dimensions
    const static int NV = 6; // uv coords of triangle corners: u0, u1, u2, v0, v1, v2
    const static int NC = 8; // C00, C01, C02, C10, C11, C12, area, anisotropy

    using VecI    = Eigen::Matrix<size_t, NV, 1>;
    using VecV    = Eigen::Matrix<double, NV, 1>;
    using VecC    = Eigen::Matrix<double, NC, 1>;
    using Triplet = Eigen::Triplet<double>;

    using MtxMap32d = Eigen::Map<Eigen::Matrix<double, 3, 2>>;
    using MtxMap23d = Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>>;

    inline double eval_f(const VecV& _x, const VecC& _c) const {

      MtxMap23d C((double*)_c.data());

      // i) Eigen way, timing is *slightly* better than the ii) version below (not much really)
      MtxMap32d X((double*)_x.data());

      //// construct S explicitly
      // Eigen::Matrix2d S;
      // S << _c[7], 0.0, 0.0, _c[7];
      // return 0.5 * _c[6] * (C * X - S).squaredNorm();

      //// construct a temp S
      return 0.5 * _c[6] * (C * X - (Eigen::Matrix2d() << _c[7], 0.0, 0.0, _c[7]).finished()).squaredNorm();

      //// ii) the other way: compute by cols
      // double f(0);
      // Eigen::Vector2d v = (C * _x.segment(0, 3) - Eigen::Vector2d(_c[7], 0));
      // f += v.squaredNorm();
      // v = C * _x.segment(3, 3) - Eigen::Vector2d(0, _c[7]);
      // f += v.squaredNorm();
      // return 0.5 * _c[6] * f;
    }

    inline void eval_gradient(const VecV& _x, const VecC& _c, VecV& _g) const {
      MtxMap23d C((double*)_c.data());
      Eigen::Matrix3d CtC = C.transpose() * C;
      _g.segment(0, 3)    = _c[6] * (CtC * _x.segment(0, 3) - C.transpose() * Eigen::Vector2d(_c[7], 0));
      _g.segment(3, 3)    = _c[6] * (CtC * _x.segment(3, 3) - C.transpose() * Eigen::Vector2d(0, _c[7]));
    }

    inline void eval_hessian(const VecV& _x, const VecC& _c, std::vector<Triplet>& _triplets) const {
      _triplets.clear();

      MtxMap23d C((double*)_c.data());
      Eigen::Matrix3d CtC = C.transpose() * C;
      CtC.array() *= _c[6];

      _triplets.clear();
      for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j) {
          _triplets.push_back(Triplet(i, j, CtC(i, j)));
          _triplets.push_back(Triplet(i + 3, j + 3, CtC(i, j)));
        }
    }

    inline double max_feasible_step(const VecV& _x, const VecV& _v, const VecC& _c) const { return DBL_MAX; }

    // generate required constants for a given input triangle
    static void constants_from_triangle( //
        const Eigen::Vector3d& _p0,      // first vertex position
        const Eigen::Vector3d& _p1,      // second vertex position
        const Eigen::Vector3d& _p2,      // third vertex position
        const Eigen::Vector3d& _fu,      // first frame vector
        const Eigen::Vector3d& _fv,      // second frame vector
        double _fscale_u,                // face scale u
        double _fscale_v,                // face scale v
        int _flabel,                     // label for this face
        double _alpha,                   // anisotropy alpha
        VecC& _c) {

      // get square-root of alpha
      const double alpha_sqrt = std::sqrt(_alpha);

      // get edge vectors
      const Eigen::Vector3d e0 = _p2 - _p1;
      const Eigen::Vector3d e1 = _p0 - _p2;
      const Eigen::Vector3d e2 = _p1 - _p0;

      // get normal vector
      Eigen::Vector3d n = -e2.cross(e1);

      // get twice the area
      double A2 = sqrt(n.dot(n));

      // normalize
      n /= A2;

      if (!std::isfinite(n.dot(n))) {
        std::cerr << "Warning: PoissonElement::constants_from_triangle had numerical issues \n";
        _c << 0, 0, 0, 0, 0, 0;
        return;
      }

      // gradient
      Eigen::Matrix3d G;
      G.col(0) = n.cross(e0);
      G.col(1) = n.cross(e1);
      G.col(2) = n.cross(e2);
      G.array() /= A2; // normalize by (2 * area)

      // frame matrix
      Eigen::Matrix<double, 3, 2> F;
      F.col(0) = _fu / _fscale_u;
      F.col(1) = _fv / _fscale_v;

      // compute and store Ft * G
      Eigen::Matrix<double, 2, 3> C = F.transpose() * G;
      if (_flabel == 2)
        C.row(0) *= alpha_sqrt;
      else if (_flabel == 1)
        C.row(1) *= alpha_sqrt;

      for (unsigned int row = 0; row < 2; ++row) // ignore last row (it's all 0)
        for (unsigned int col = 0; col < 3; ++col)
          _c[3 * row + col] = C(row, col); // order is 0: C00, 1: C01, 2: C02 (row 0), 3: C10, 4: C11, 5: C12 (row 1)

      _c[6] = 0.5 * A2;
      _c[7] = alpha_sqrt;
    }
  };

  //===========================================================================
  //== CLASS DEFINITION =======================================================
  //===========================================================================

  class SnappingElement {
  public:
    const static int NV = 2; // 0: real coordinate, 1: auxiliary integer variable
    const static int NC = 1; // snapping weight

    using VecI    = Eigen::Matrix<size_t, NV, 1>;
    using VecV    = Eigen::Matrix<double, NV, 1>;
    using VecC    = Eigen::Matrix<double, NC, 1>;
    using Triplet = Eigen::Triplet<double>;

    inline double eval_f(const VecV& _x, const VecC& _c) const { //
      return 0.5 * _c[0] * (_x[0] - _x[1]) * (_x[0] - _x[1]);
    }

    inline void eval_gradient(const VecV& _x, const VecC& _c, VecV& _g) const {
      _g[0] = _c[0] * (_x[0] - _x[1]);
      _g[1] = _c[0] * (_x[1] - _x[0]);
    }

    inline void eval_hessian(const VecV& _x, const VecC& _c, std::vector<Triplet>& _triplets) const {
      _triplets.clear();
      _triplets.push_back(Triplet(0, 0, _c[0]));
      _triplets.push_back(Triplet(1, 1, _c[0]));
      _triplets.push_back(Triplet(0, 1, -_c[0]));
      _triplets.push_back(Triplet(1, 0, -_c[0]));
    }

    inline double max_feasible_step(const VecV& _x, const VecV& _v, const VecC& _c) const { return DBL_MAX; }
  };

  // ===========================================================================

  class RegularizationElement {

  public:
    // define dimensions
    const static int NV = 1; // uv coordinate
    const static int NC = 1; // reg weight

    using VecI    = Eigen::Matrix<size_t, NV, 1>;
    using VecV    = Eigen::Matrix<double, NV, 1>;
    using VecC    = Eigen::Matrix<double, NC, 1>;
    using Triplet = Eigen::Triplet<double>;

    inline double eval_f(const VecV& _x, const VecC& _c) const { return 0.5 * _c[0] * _x[0] * _x[0]; }

    inline void eval_gradient(const VecV& _x, const VecC& _c, VecV& _g) const { _g[0] = _c[0] * _x[0]; }

    inline void eval_hessian(const VecV& _x, const VecC& _c, std::vector<Triplet>& _triplets) const {
      _triplets.clear();
      _triplets.push_back(Triplet(0, 0, _c[0]));
    }

    inline double max_feasible_step(const VecV& _x, const VecV& _v, const VecC& _c) const { return DBL_MAX; }
  };

} // namespace IGSV

#pragma once
