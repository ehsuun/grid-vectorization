#pragma once

#include "_cgal_delaunay_2.h" // using IndexedDT
#include <Eigen/Core>
#include <opencv2/opencv.hpp>

namespace IGSV {

  class FaceData;

  bool trace_streamlines(const Eigen::MatrixXi& _TT,                                //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isBlack,  //
                         const Eigen::MatrixXd& _F_barycenters,                     //
                         const Eigen::VectorXcd (&_F_frameField)[4],                //
                         const Eigen::MatrixXi& _E_periodJumps,                     //
                         const IndexedDT& _tri,                                     //
                         const cv::Mat& _bw,                                        //
                         double _sw_avg,                                            //
                         double _scale_multiplier,                                  //
                         double _tangent_ratio_threshold,                           //
                         double _tangent_ratio_streamlen,                           //
                         Eigen::VectorXi& _F_labels,                                //
                         std::vector<FaceData>& _FaceDataList);

}; // namespace IGSV
