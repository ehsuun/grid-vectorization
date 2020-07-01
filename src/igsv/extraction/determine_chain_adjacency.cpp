#include <igsv/extraction/determine_chain_adjacency.h>

#include <igsv/extraction_entities/Chain.h>
#include <igsv/extraction_entities/QEdge.h>

namespace IGSV {

  bool determine_chain_adjacency(const std::vector<QEdge>& _QE,                        //
                                 const std::vector<Chain>& _CH,                        //
                                 std::vector<std::pair<int, int>>& _CH_CH_connections, //
                                 std::vector<bool>& _CH_CH_connections_enabled) {

    Eigen::MatrixXi A_chains;
    A_chains.setZero(_CH.size(), _CH.size());

    int c0, c1;
    for (int chi0 = 0; chi0 < _CH.size(); ++chi0) {

      // skip if this chain was removed
      if (_CH[chi0].is_removed)
        continue;

      for (auto qei : _CH[chi0].adjacent_qedges) {

        if (!_QE[qei].is_used)
          continue;

        const int chi1 = _QE[qei].chain;
        if (chi1 == -1)
          continue;

        // lower index always first
        c0 = chi0 < chi1 ? chi0 : chi1;
        c1 = chi0 < chi1 ? chi1 : chi0;

        // skip if the connection already exists
        if (A_chains(c0, c1) == 1)
          continue;

        // otherwise, connect
        A_chains(c0, c1) = 1;
        A_chains(c1, c0) = 1;
        _CH_CH_connections.emplace_back(c0, c1);
      }
    }

    std::sort(_CH_CH_connections.begin(), _CH_CH_connections.end());
    _CH_CH_connections.erase(std::unique(_CH_CH_connections.begin(), _CH_CH_connections.end()),
                             _CH_CH_connections.end());

    _CH_CH_connections_enabled.resize(_CH_CH_connections.size());
    std::fill(_CH_CH_connections_enabled.begin(), _CH_CH_connections_enabled.end(), false);

    return true;
  }

} // namespace IGSV
