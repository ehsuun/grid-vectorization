#include <igsv/extraction_entities/Bezier.h>

namespace IGSV {
  namespace Bezier {

    //===========================================================================

    double Bernstein(int n, int k, double t) { //
      return Pascal[n][k] * std::pow(t, k) * std::pow(1 - t, n - k);
    }

    //===========================================================================

    std::complex<double> Casteljau(const Eigen::Vector4cd& CP, int k, int i, double t) {
      return (k == 0) ? CP(i)
                      : (1.0 - t) * Casteljau(CP, k - 1, i, t) //
                            + t * Casteljau(CP, k - 1, i + 1, t);
    }

    //===========================================================================

    std::complex<double> BezierPoint(const Eigen::Vector4cd& CP, double t) { //
      return Casteljau(CP, 3, 0, t);
    }

    //===========================================================================

  } // namespace Bezier
} // namespace IGSV
