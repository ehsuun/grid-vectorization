#pragma once

namespace IGSV {

  bool compute_scale(const Eigen::MatrixXd& _BC,  //
                     const cv::Mat& _detail_mask, //
                     double _mask_factor,         //
                     double _avg_stroke_width,    //
                     double _scale_multiplier,    //
                     Eigen::VectorXd& _F_scale);

}
