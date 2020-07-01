//
// Created by Tibor Stanko on 10/04/2020.
//

#ifndef IGSV_COLORS_H
#define IGSV_COLORS_H

#include <opencv2/opencv.hpp>
namespace IGSV {
  namespace colors {
    const cv::Scalar orange(64, 163, 241, 255);
    const cv::Scalar violet(195, 142, 153, 255);
    const cv::Scalar white(255, 255, 255, 255);
    const cv::Scalar black(0, 0, 0, 255);
    const cv::Scalar transparent(255, 255, 255, 0);
  } // namespace colors
} // namespace IGSV

#endif // IGSV_COLORS_H
