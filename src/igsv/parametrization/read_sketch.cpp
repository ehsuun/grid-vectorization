#include <igsv/parametrization/read_sketch.h>

#include <igsv/common/defs.h>
#include <igsv/common/print_timings.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace IGSV {

  // ===========================================================================

  cv::Mat strel(const long k) {
    return cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(2 * k + 1, 2 * k + 1), cv::Point(k, k));
  }

  // ===========================================================================

  bool read_sketch(const std::string& _filename_input, //
                   const std::string& _filename_mask,  //
                   float& _mask_factor,                //
                   cv::Mat& _input,                    //
                   cv::Mat& _detail_mask,              //
                   bool& _mask_read_from_file) {

    //// read input sketch
    _input = cv::imread(_filename_input, cv::IMREAD_COLOR);

    if (_input.cols == 0 || _input.rows == 0) {
      // failed to read the sketch! abort.
      log_error(std::string("could not read the input from ") + _filename_input);
      return false;
    }

    // print info
    log_message((boost::format("input read from %s [%dx%d]") % _filename_input % _input.cols % _input.rows).str());

    //// read the mask
    _mask_read_from_file = false;

    if (!_filename_mask.empty()) {
      if (_mask_factor == 0.0f) {
        log_warning("mask factor set to 0, ignoring the specified mask [" + _filename_mask + "]");

      } else {
        // try to read the mask
        _detail_mask = cv::imread(_filename_mask, cv::IMREAD_GRAYSCALE);

        if (_detail_mask.rows == 0 || _detail_mask.cols == 0) {
          // mask not read, continue without it
          log_warning(std::string("mask not read from ") + _filename_mask);

        } else {
          // mask read ok
          if (_input.rows != _detail_mask.rows || _input.cols != _detail_mask.cols) {

            // invalid mask
            log_warning((boost::format("mask resolution [%dx%d] does not match input resolution [%dx%d], ignoring") //
                         % _detail_mask.cols                                                                        //
                         % _detail_mask.rows                                                                        //
                         % _input.cols                                                                              //
                         % _input.rows)
                            .str());

          } else {
            // mask is valid
            _mask_read_from_file = true;

            // print info
            log_message(
                (boost::format("mask read from %s [%dx%d]") % _filename_mask % _detail_mask.cols % _detail_mask.rows)
                    .str());
          }
        }
      }
    }

    // reset the mask if not read (set to grey)
    if (!_mask_read_from_file) {
      _detail_mask = cv::Mat::ones(_input.size(), CV_8UC1);
      _detail_mask.setTo(128);
    }

    return true;
  }

  // ===========================================================================

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
                      double& _sw_stddev) {

    // 1. Convert to greyscale
    cv::cvtColor(_input, _grey, cv::COLOR_BGR2GRAY);
    _grey.convertTo(_grey, CV_8UC1);

    // 2. Binarize
    cv::threshold(_grey, _bw, _threshold, 255, cv::THRESH_BINARY);
    cv::normalize(_bw, _bw, 0, 255, cv::NORM_MINMAX, CV_8UC1);

    // 3. Distance transform on the binarized sketch
    cv::distanceTransform(_bw, _dt, cv::DIST_L1, cv::DIST_MASK_5);
    _dt.convertTo(_dt, CV_8UC1);

    // 4. Distance transform on the inverted binarized sketch
    cv::Mat idt;
    cv::distanceTransform(255 - _bw, idt, cv::DIST_L1, cv::DIST_MASK_5);
    idt.convertTo(idt, CV_8UC1);
    cv::normalize(idt, _idt_01, 0.0, 1.0, cv::NORM_MINMAX, CV_64F); // normalized

    // 5. Stroke width
    double idt_min, idt_max;
    cv::minMaxLoc(idt, &idt_min, &idt_max);

    _sw = idt.clone();
    cv::dilate(_sw, _sw, strel(std::round<long>(idt_max)));
    _sw.convertTo(_sw, CV_64F);
    _sw *= 2; // multiply by a factor of 2 to get full width from the distance?

    double sw_min;
    cv::minMaxLoc(_sw, &sw_min, &_sw_max);
    cv::normalize(_sw, _sw_01, 0.0, 1.0, cv::NORM_MINMAX, CV_64F); // normalized

    // 6. Stroke width stats
    cv::Mat sw_mask;
    cv::threshold(_sw, sw_mask, 0, 1, cv::THRESH_BINARY);
    sw_mask.convertTo(sw_mask, CV_8UC1);
    cv::Scalar sw_avg, sw_stddev;
    cv::meanStdDev(_sw, sw_avg, sw_stddev, sw_mask);
    _sw_avg    = sw_avg[0];
    _sw_stddev = sw_stddev[0];

    DEBUG_PRINT_INFO("stroke width stats | min=%0.2f | max=%0.2f | avg=%0.2f | stddev=%0.2f |", //
                     sw_min, _sw_max, _sw_avg, _sw_stddev);

    cv::Mat sw_trim = _sw.clone();
    sw_trim.convertTo(sw_trim, CV_32F);
    sw_trim.setTo(_sw_avg + 3 * _sw_stddev, _sw == 0);
    cv::GaussianBlur(sw_trim, _sw, cv::Size(21, 21), _sw_avg);

    return true;
  }

  // ===========================================================================

} // namespace IGSV
