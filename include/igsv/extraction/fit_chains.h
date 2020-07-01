//=============================================================================
//
//  FUNCTION : fit_chains
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <map>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  struct QVertex;
  struct QEdge;
  struct Chain;

  //== FUNCTION DEFINITION ====================================================

  bool fit_chains(const std::vector<QVertex>& _QV,                         //
                  const std::map<std::pair<int, int>, int>& _QVPair_to_QE, //
                  std::vector<QEdge>& _QE,                                 //
                  std::vector<Chain>& _CH,                                 //
                  double _avg_stroke_width,                                //
                  double _narrow_band_radius);

} // namespace IGSV
