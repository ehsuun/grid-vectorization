//=============================================================================
//
//  FUNCTION : assign_samples_to_chains
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  struct QEdge;
  struct Pixel;
  struct Chain;

  //== FUNCTION DEFINITION ====================================================

  bool assign_samples_to_chains(const std::vector<QEdge>& _QE, //
                                std::vector<Pixel>& _PX,       //
                                std::vector<Chain>& _CH);

} // namespace IGSV
