#include <igsv/extraction/init_base_labels.h>
#include <igsv/common/defs.h>
#include <igsv/extraction_entities/QEdge.h>
#include <igsv/extraction_entities/QPort.h>
#include <igsv/extraction_entities/QVertex.h>
#include <igsv/extraction/connected_components.h>

#include <iostream>

namespace IGSV {

  bool init_base_labels(const std::vector<std::pair<int, int>>& _QP_QP_connections, //
                        const bool split_at_nodes,                                  //
                        std::vector<QVertex>& _QV,                                  //
                        std::vector<QPort>& _QP,                                    //
                        std::vector<QEdge>& _QE) {

    const int n_qverts = _QV.size();
    const int n_qports = _QP.size();
    const int n_qedges = _QE.size();

    // Two q-edges are in the same chain if locally, they were extracted from the same isoline.
    std::vector<std::vector<int>> QP_config;
    if (split_at_nodes) {
      std::vector<std::pair<int, int>> QP_QP_active;
      int qpi0, qpi1, qvi0, qvi1;
      for (auto qpi_pair : _QP_QP_connections) {
        qpi0 = qpi_pair.first;
        qpi1 = qpi_pair.second;
        qvi0 = _QP[qpi0].qvi;
        qvi1 = _QP[qpi1].qvi;
        if (_QP[qpi0].is_irregular || _QV[qvi0].is_split_node || //
            _QP[qpi1].is_irregular || _QV[qvi1].is_split_node)
          continue;
        if (_QV[qvi0].valence < 3 || _QV[qvi1].valence < 3)
          QP_QP_active.push_back(qpi_pair);
      }
      connected_components(QP_QP_active, QP_config, n_qports);
    } else {
      connected_components(_QP_QP_connections, QP_config, n_qports);
    }

    int n_chains = 0;
    int qei;

    for (int cc = 0; cc < QP_config.size(); cc++) {
      bool use_this_component = false;

      // is at least one edge in the chain tangent?
      for (int qpi : QP_config[cc]) {
        qei = _QP[qpi].qei;

        // skip if the q-port is not connected to an edge
        if (qei == -1)
          continue;

        if (split_at_nodes) {
          // skip if the q-edge is not used
          if (!_QE[qei].is_used)
            continue;
        } else {
          // skip if the q-edge is not tangent
          if (!_QE[qei].is_tangent)
            continue;
        }

        // found a valid tangent q-edge, exit the loop
        use_this_component = true;
        break;
      }

      if (!use_this_component)
        continue; // move on to the next component if this one does not have a tangent edge

      for (int qpi : QP_config[cc]) {
        _QP[qpi].cc = cc;

        qei = _QP[qpi].qei;

        if (qei == -1)
          continue;

        _QE[qei].base_label = n_chains;
        if (!split_at_nodes)
          _QE[qei].is_used = true;
      }
      n_chains++;
    }

    DEBUG_PRINT_INFO("%d chains", n_chains);

    if (!split_at_nodes) {
      for (int qvi = 0; qvi < n_qverts; ++qvi)
        _QV[qvi].valence = 0;

      for (int qei = 0; qei < n_qedges; ++qei)
        if (_QE[qei].is_used) {
          _QV[_QE[qei].qvi0].valence++;
          _QV[_QE[qei].qvi1].valence++;
        }
    }

    return true;
  }

} // namespace IGSV
