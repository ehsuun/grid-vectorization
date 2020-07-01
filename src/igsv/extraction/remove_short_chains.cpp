#include <igsv/extraction/remove_short_chains.h>

#include <igsv/common/defs.h>
#include <igsv/extraction_entities/Chain.h>
#include <igsv/extraction_entities/QEdge.h>
#include <igsv/extraction_entities/QVertex.h>

#define HISTOGRAM_SIZE 4

namespace IGSV {

  bool remove_short_chains(std::vector<QVertex>& _QV, //
                           std::vector<QEdge>& _QE,   //
                           std::vector<Chain>& _CH) {

    int n_to_remove = 0;
    int chi         = 0;

    for (auto& _chain : _CH) {

      //// Determine the parametric coverage of this chain
      //// via a histogram (# bins = HISTOGRAM_SIZE)
      _chain.histogram.resize(HISTOGRAM_SIZE, 1);
      _chain.histogram.fill(0);
      for (auto const& _sample : _chain.samples) {
        for (int k = 1; k <= HISTOGRAM_SIZE; k++) {
          if (k == HISTOGRAM_SIZE) {
            // we have arrived at the end : increase count in the last bin
            _chain.histogram[k - 1]++;
            break;
          } else {
            // check the next bin : increase the count if t_chain is in it (meaning it wasn't in the previous bins)
            if (_sample.t_chain < (double)k / (double)HISTOGRAM_SIZE) {
              _chain.histogram[k - 1]++;
              break;
            }
          }
        }
      }

      //// If any one bin is empty, we don't have enough samples for fitting -> discard the chain
      if ((_chain.histogram.array() == 0).any()) {
        DEBUG_PRINT_WARNING("chain #%d has too few samples assigned, remove", chi);

        // for (int k = 0; k < HISTOGRAM_SIZE; k++) {
        //   DEBUG_PRINT_INFO("    [%0.2f, %0.2f] : %d samples %s",     //
        //                    (double)k / (double)HISTOGRAM_SIZE,       //
        //                    (double)(k + 1) / (double)HISTOGRAM_SIZE, //
        //                    _chain.histogram[k],                             //
        //                    _chain.histogram[k] == 0 ? " (!)" : "");
        // }

        if (_chain.qedges.size() == 1) {
          int qei = *(_chain.qedges.begin());

          if (_QE[qei].has_flipped_faces) {
            DEBUG_PRINT_WARNING("chain #%d, cannot disable edge %d, which has flipped faces", chi, qei);

          } else {
            DEBUG_PRINT_WARNING("chain #%d, disable edge %d", chi, qei);
            _QE[qei].base_label = -1;
            _QE[qei].is_used    = false;
            _QE[qei].chain      = -1;

            _chain.is_removed = true;
            n_to_remove++;
          }

        } else {
          const int middle_vertex          = _chain.qverts[_chain.qverts.size() / 2];
          _QV[middle_vertex].is_split_node = true;
          DEBUG_PRINT_WARNING("chain #%d, split in the middle (qv%d)", chi, middle_vertex);
        }
      }

      chi++;
    }

    DEBUG_PRINT_INFO("removing %d chains", n_to_remove);

    return true;
  }

} // namespace IGSV
