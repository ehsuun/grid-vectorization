//=============================================================================
//
//  STRUCT : Pixel
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <vector>

#include "Sample.h"

//== NAMESPACES ===============================================================

namespace IGSV {

  //== STRUCT DEFINITION ======================================================

  struct Pixel {
    double f0x, f0y, f1x, f1y; // (combed) frame at the pixel

    std::vector<Sample> samples; // list of samples, ordered by distance
    bool has_active_samples;
  };

} // namespace IGSV
