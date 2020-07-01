#include <igsv/extraction/assign_samples_to_chains.h>

#include <igsv/extraction_entities/Chain.h>
#include <igsv/extraction_entities/Pixel.h>
#include <igsv/extraction_entities/QEdge.h>

#include <igsv/common/defs.h>

namespace IGSV {

  bool assign_samples_to_chains(const std::vector<QEdge>& _QE, //
                                std::vector<Pixel>& _PX,       //
                                std::vector<Chain>& _CH) {

    int n_pixels_without_active_samples = 0;

    // temp vars
    int pid, qei, chi, qvi0, qvi1, i0, i1;
    double t_chain;

    pid = 0;

    for (auto& pixel : _PX) {

      pixel.has_active_samples = false;

      // loop over samples for this pixel
      for (auto const& sample : pixel.samples) {

        qei = sample.qei;

        // skip if the edge is not used
        if (!_QE[qei].is_used)
          continue;

        chi = _QE[qei].chain;

        // skip if no chain is assigned to the edge
        if (chi == -1)
          continue;

        // skip if this chain was removed
        if (_CH[chi].is_removed)
          continue;

        // global vertex indices
        qvi0 = _QE[qei].qvi0;
        qvi1 = _QE[qei].qvi1;

        // local vertex indices w.r.t. the chain
        i0 = _CH[chi].vindex_map[qvi0];
        i1 = _CH[chi].vindex_map[qvi1];

        // parameter w.r.t. the chain
        t_chain = (1.0 - sample.t) * _CH[chi].t[i0] + sample.t * _CH[chi].t[i1];

        // assign this pixel sample to the chain
        _CH[chi].samples.emplace_back(sample, chi, t_chain);

        // ignore further samples for this pixel
        pixel.has_active_samples = true;
        break;

      } // loop over pixel samples

      if (!pixel.has_active_samples)
        n_pixels_without_active_samples++;

      pid++;
    }

    if (n_pixels_without_active_samples > 0)
      DEBUG_PRINT_WARNING("%d pixels without active samples", n_pixels_without_active_samples);

    return true;
  }
} // namespace IGSV
