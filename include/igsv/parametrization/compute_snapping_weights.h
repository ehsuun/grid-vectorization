#pragma once

namespace IGSV {

  bool compute_snapping_weights(const Eigen::MatrixXd& _V, //
                                const cv::Mat& _idt_01,    //
                                Eigen::VectorXd& _V_snappingWeights);

}
