//=============================================================================
//
//  FUNCTION : assign_pixels_to_qedges
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <opencv2/opencv.hpp>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  struct QVertex;
  struct QPort;
  struct QEdge;
  struct Pixel;

  //== FUNCTION DEFINITION ====================================================

  bool assign_pixels_to_qedges(const cv::Mat& _grey,                 //
                               const cv::Mat& _sw_01,                //
                               const Eigen::VectorXcd& _PX_XY,       //
                               const Eigen::VectorXi& _PX_I,         //
                               const Eigen::VectorXi& _PX_J,         //
                               const Eigen::VectorXi& _PX_fid,       //
                               const Eigen::MatrixXd& _X0_unit_comb, //
                               const Eigen::MatrixXd& _X1_unit_comb, //
                               const std::vector<QVertex>& _QV,      //
                               const std::vector<QPort>& _QP,        //
                               const std::vector<QEdge>& _QE,        //
                               const double _sw_avg,                 //
                               std::vector<Pixel>& _PX);

} // namespace IGSV
