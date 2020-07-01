//=============================================================================
//
//  FUNCTION : parametrize_chains
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  struct QVertex;
  struct QEdge;
  struct Chain;

  //== FUNCTION DEFINITION ====================================================

  bool parametrize_chains(const std::vector<QVertex>& _QV, //
                          const std::vector<QEdge>& _QE,   //
                          std::vector<Chain>& _CH);

} // namespace IGSV
