//=============================================================================
//
//  FUNCTION : determine_chain_adjacency
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  struct QEdge;
  struct Chain;

  //== FUNCTION DEFINITION ====================================================

  bool determine_chain_adjacency(const std::vector<QEdge>& _QE,                        //
                                 const std::vector<Chain>& _CH,                        //
                                 std::vector<std::pair<int, int>>& _CH_CH_connections, //
                                 std::vector<bool>& _CH_CH_connections_enabled);

} // namespace IGSV
