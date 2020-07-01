//=============================================================================
//
//  CLASS : FaceData
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <vector>

//== FORWARD DECLARATIONS =====================================================

//== NAMESPACES ===============================================================

namespace IGSV {

  //== TYPEDEFS ===============================================================

  //== CLASS DEFINITION =======================================================

  class FaceData {
  public:
    FaceData(int _f, double _bx, double _by, int _numsamples = 0) //
        : f(_f),                                                  //
          bx(_bx),                                                //
          by(_by),                                                //
          numsamples(_numsamples) {

      for (int k = 0; k < 4; k++) {
        samples[k].resize(0, 3);
        fpath[k].clear();
        epath[k].clear();
      }

      c[0] = 0;
      c[1] = 0;
    }

    // face index
    int f;

    // barycenter position
    double bx, by;

    // black pixel counts (two directions)
    int c[2];

    // number of samples in four directions (excluding the barycenter)
    int numsamples;

    // traced streamlines
    Eigen::MatrixXd samples[4];

    // triangles and edges intersected by traced streamlines
    std::vector<int> fpath[4];
    std::vector<int> epath[4];
  };
} // namespace IGSV
