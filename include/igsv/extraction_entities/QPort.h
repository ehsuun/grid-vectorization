//=============================================================================
//
//  STRUCT : QPort
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== STRUCT DEFINITION ======================================================
  struct QPort {

    QPort(int _qvi, int _pj, int _f, int _u, int _v, double _x_end, double _y_end)
        : qvi(_qvi), pj(_pj), f(_f), u(_u), v(_v), x_end(_x_end), y_end(_y_end) {}

    // const properties
    const int qvi;      // index of the associated qvertex
    const int pj;       // period jump
    const int f;        // index of the face the qport point into
    const int u;        // u-coords of the corresponding corner
    const int v;        // v-coord of the corresponding corner
    const double x_end; // x-coord of the endpoint, for visualisation
    const double y_end; // y-coord of the endpoint, for visualisation

    // modifiable properties
    int qei           = -1;    // index of the associated qedge (-1 if no edge
    int cc            = -1;    // index of the connected component
    bool is_irregular = false; // is the qport irregular?
  };

} // namespace IGSV
