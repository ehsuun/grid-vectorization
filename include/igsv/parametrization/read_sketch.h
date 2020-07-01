#pragma once

#include <opencv2/opencv.hpp>

namespace IGSV {

  bool read_sketch(const std::string& _filename_input, //
                   const std::string& _filename_mask,  //
                   float& _mask_factor,                //
                   cv::Mat& _input,                    //
                   cv::Mat& _detail_mask,              //
                   bool& _mask_read_from_file);

  bool process_sketch(const cv::Mat& _input,       //
                      const cv::Mat& _detail_mask, //
                      int _threshold,              //
                      cv::Mat& _grey,              //
                      cv::Mat& _bw,                //
                      cv::Mat& _dt,                //
                      cv::Mat& _idt_01,            //
                      cv::Mat& _sw,                //
                      cv::Mat& _sw_01,             //
                      double& _sw_max,             //
                      double& _sw_avg,             //
                      double& _sw_stddev);

} // namespace IGSV
