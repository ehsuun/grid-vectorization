//=============================================================================
//
//  FUNCTION : compute_narrow_band
//
//=============================================================================

#pragma once

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================

  bool compute_narrow_band(const Eigen::MatrixXd& _V,                                   //
                           const Eigen::MatrixXi& _F,                                   //
                           const Eigen::MatrixXi& _TT,                                  //
                           const Eigen::MatrixXi& _TTi,                                 //
                           const std::vector<std::vector<int>>& _VT,                    //
                           const std::vector<std::vector<int>>& _VTi,                   //
                           const Eigen::MatrixXd& _BC,                                  //
                           const cv::Mat& _bw,                                          //
                           const cv::Mat& _dt,                                          //
                           const cv::Mat& _detail_mask,                                 //
                           double _mask_factor,                                   //
                           double _narrow_band_radius,                            //
                           double _avg_stroke_width,                              //
                           Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isBlack,          //
                           Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrow,         //
                           Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrowBoundary, //
                           Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow,         //
                           Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrow,         //
                           Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrowBoundary);

} // namespace IGSV
