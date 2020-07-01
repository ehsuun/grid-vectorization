//=============================================================================
//
//  FUNCTION : connected_components
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================
  // based on https://stackoverflow.com/questions/46659904/connected-components-in-an-undirected-graph-in-c

  void connected_components(const Eigen::MatrixXi& E, std::vector<std::vector<int>>& conncomp, int n);
  void connected_components(const std::vector<std::pair<int, int>>& edges, std::vector<std::vector<int>>& conncomp,
                            int n);

} // namespace IGSV
