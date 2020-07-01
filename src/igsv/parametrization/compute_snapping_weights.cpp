#include <Eigen/Core>
#include <opencv2/opencv.hpp>

#include <igsv/parametrization/compute_snapping_weights.h>
#include <igsv/parametrization/pix_coords_from_xy.h>

namespace IGSV {

  bool compute_snapping_weights(const Eigen::MatrixXd& _V, //
                                const cv::Mat& _idt_01,    //
                                Eigen::VectorXd& _V_snappingWeights) {

    const double sigma    = 0.1;
    const double sigma_sq = sigma * sigma;

    int i, j;

    _V_snappingWeights.resize(_V.rows(), 1);
    _V_snappingWeights.fill(1.0);

    double weight;
    for (int v = 0; v < _V.rows(); ++v) {

      // get pixel coordinates for this vertex
      pix_coords_from_xy(_idt_01.cols, _idt_01.rows, _V(v, 0), _V(v, 1), i, j);

      // inverse distance weight (gaussian, max=1.0 on the stroke centerline)
      weight = _idt_01.at<double>(i, j);
      weight = std::exp(-0.5 * weight * weight / sigma_sq);
      if (weight < 0.)
        weight = 0.;
      if (weight > 1.)
        weight = 1.;
      _V_snappingWeights(v) *= (1.0 - weight);
    }

    // clamp to [0, 1]
    _V_snappingWeights = _V_snappingWeights.cwiseMax(0.0);
    _V_snappingWeights = _V_snappingWeights.cwiseMin(1.0);

    return true;
  }
} // namespace IGSV
