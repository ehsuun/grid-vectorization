//=============================================================================
//
//  FUNCTION : init_base_labels
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  struct QVertex;
  struct QPort;
  struct QEdge;

  //== FUNCTION DEFINITION ====================================================

  bool init_base_labels(const std::vector<std::pair<int, int>>& _QP_QP_connections, //
                        const bool split_at_nodes,                                  //
                        std::vector<QVertex>& _QV,                                  //
                        std::vector<QPort>& _QP,                                    //
                        std::vector<QEdge>& _QE);

} // namespace IGSV
