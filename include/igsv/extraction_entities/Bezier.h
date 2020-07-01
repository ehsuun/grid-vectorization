//=============================================================================
//
//  FUNCTION : Bernstein
//             Casteljau
//             BezierPoint
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <cmath>
#include <complex>

//== NAMESPACES ===============================================================

namespace IGSV {
  namespace Bezier {

    //===========================================================================

    const int Pascal[4][4] = {
      { 1, 0, 0, 0 }, //
      { 1, 1, 0, 0 }, //
      { 1, 2, 1, 0 }, //
      { 1, 3, 3, 1 }  //
    };

    //===========================================================================

    double Bernstein(int n, int k, double t);
    std::complex<double> Casteljau(const Eigen::Vector4cd& CP, int k, int i, double t);
    std::complex<double> BezierPoint(const Eigen::Vector4cd& CP, double t);

  } // namespace Bezier
} // namespace IGSV
