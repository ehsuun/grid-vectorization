//=============================================================================
//
//  FUNCTION : init_chains
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

  bool init_chains(const std::vector<QVertex>& _QV, //
                   std::vector<QEdge>& _QE,         //
                   std::vector<Chain>& _CH);

} // namespace IGSV
