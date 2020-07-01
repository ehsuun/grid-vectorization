#include <igsv/parametrization/pix_coords_from_xy.h>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>

namespace IGSV {

  bool compute_scale(const Eigen::MatrixXd& _BC,  //
                     const cv::Mat& _detail_mask, //
                     double _mask_factor,         //
                     double _avg_stroke_width,    //
                     double _scale_multiplier,    //
                     Eigen::VectorXd& _F_scale) {

    cv::Mat _scale_mask = cv::Mat::ones(_detail_mask.size(), CV_32F);
    if (_mask_factor > 1.0f) {
      // white pixels : coarser scale
      _scale_mask.setTo(_mask_factor, _detail_mask == 255);

      // black pixels : finer scale
      _scale_mask.setTo(1.0f / _mask_factor, _detail_mask == 0);
    }

    int i, j;

    _F_scale.setZero(_BC.rows(), 1);
    for (int f = 0; f < _BC.rows(); ++f) {
      // get pixel coordinates of triangle's barycenter
      pix_coords_from_xy(_scale_mask.cols, _scale_mask.rows, _BC(f, 0), _BC(f, 1), i, j);

      // divide by scale mask
      _F_scale(f) = 1.0 / _scale_mask.at<float>(i, j);
    }

    // normalize by the average stroke width and the global scale multiplier
    _F_scale.array() /= (_avg_stroke_width * _scale_multiplier);

    return true;
  }

} // namespace IGSV
