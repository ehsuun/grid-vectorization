#pragma once

#include <Eigen/Core>
#include <opencv2/opencv.hpp>

namespace IGSV {

  bool triangulate(const cv::Mat& _bw,                  //
                   const cv::Mat& _dt,                  //
                   double _maxarea,                     //
                   bool _adaptive,                      //
                   double _narrow_band_radius,          //
                   double _avg_stroke_width,            //
                   Eigen::MatrixXd& _V,                 //
                   Eigen::MatrixXi& _F,                 //
                   Eigen::MatrixXd& _UV_generated,      //
                   Eigen::MatrixXd& _BC,                //
                   Eigen::MatrixXi& _TT,                //
                   Eigen::MatrixXi& _TTi,               //
                   std::vector<std::vector<int>>& _VT,  //
                   std::vector<std::vector<int>>& _VTi, //
                   Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isBlack);

  bool process_mesh(const cv::Mat& _bw,                  //
                    const Eigen::MatrixXd& _V,           //
                    const Eigen::MatrixXi& _F,           //
                    Eigen::MatrixXd& _UV_generated,      //
                    Eigen::MatrixXd& _BC,                //
                    Eigen::MatrixXi& _TT,                //
                    Eigen::MatrixXi& _TTi,               //
                    std::vector<std::vector<int>>& _VT,  //
                    std::vector<std::vector<int>>& _VTi, //
                    Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isBlack);

}; // namespace IGSV
