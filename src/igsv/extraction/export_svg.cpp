#include <igsv/common/print_timings.h>
#include <igsv/extraction_entities/Chain.h>

// includes needed by simple_svg
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <simple_svg/simple_svg_1.0.0.hpp>

namespace IGSV {

  bool export_svg(const std::string& _filename_output, //
                  const std::vector<Chain>& _chains,   //
                  int _width, int _height) {

    svg::Stroke black_stroke   = svg::Stroke(std::max(std::floor(0.004 * (double)_width), 1.0), svg::Color::Black);
    svg::Fill transparent_fill = svg::Fill();

    svg::Document svg_file(_filename_output, svg::Layout(svg::Dimensions(_width, _height), svg::Layout::BottomLeft));

    for (const auto& _chain : _chains)
      for (int i = 0; i < _chain.CP_spline.cols(); ++i) {
        svg::Path bezier_curve(transparent_fill, black_stroke);
        for (int k = 0; k < 4; k++)
          bezier_curve << svg::Point(_chain.CP_spline(k, i).real() + 0.5, //
                                     _chain.CP_spline(k, i).imag() + 0.5);
        svg_file << bezier_curve;
      }

    svg_file.save();

    log_message((boost::format("saved to %s") % _filename_output).str());

    return true;
  }

} // namespace IGSV
