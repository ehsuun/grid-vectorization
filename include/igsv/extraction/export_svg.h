#pragma once

namespace IGSV {

  class Chain;

  bool export_svg(const std::string& _filename_output, //
                  const std::vector<Chain>& _chains,   //
                  int _width, int _height);

} // namespace IGSV
