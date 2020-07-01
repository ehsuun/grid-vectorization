//=============================================================================
//
//  FUNCTION : compute_constraints
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <opencv2/opencv.hpp>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================

  bool compute_constraints(
      const cv::Mat& _grey,               // grey image
      const cv::Mat& _bw,                 // binary image
      double _sw_avg,                     // average stroke width
      cv::Mat& _grey_blurred,             // grey image, blurred
      unsigned& _n_dark,                  // number of black pixels
      Eigen::MatrixXd& _pix_weight,       // weights (all pixels)
      Eigen::MatrixXi& _pix2band,         // conversion: pixel indices (i,j) -> band indices (k), _pix2band(i,j) = k
      Eigen::VectorXi& _PX_I,             // pixel coords : i
      Eigen::VectorXi& _PX_J,             // pixel coords : j
      Eigen::VectorXcd& _PX_XY,           // pixel position (complex)
      Eigen::VectorXcd& _PX_SobelTangent, // Sobel tangent (complex)
      Eigen::VectorXd& _PX_SobelWeight    // Sobel tangent confidence (scalar)
  );

} // namespace IGSV
