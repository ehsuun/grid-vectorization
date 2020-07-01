#pragma once

#include <vector>

namespace IGSV {

  struct QVertex;
  struct QEdge;
  struct Chain;
  
  bool remove_short_chains(std::vector<QVertex>& _QV, //
                           std::vector<QEdge>& _QE,   //
                           std::vector<Chain>& _CH);

} // namespace IGSV
