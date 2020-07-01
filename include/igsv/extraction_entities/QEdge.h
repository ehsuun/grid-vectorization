//=============================================================================
//
//  STRUCT : QEdge
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== STRUCT DEFINITION ======================================================
  struct QEdge {

    QEdge(int _qpi0, int _qpi1, int _qvi0, int _qvi1, double _length, bool _has_flipped_faces, bool _is_tangent)
        : qpi0(_qpi0),
          qpi1(_qpi1),
          qvi0(_qvi0),
          qvi1(_qvi1),
          length(_length),
          has_flipped_faces(_has_flipped_faces),
          is_tangent(_is_tangent) {}

    // const properties
    const int qpi0; // first qport
    const int qpi1; // second qport
    const int qvi0; // first qvertex
    const int qvi1; // second qvertex
    const double length;
    const bool has_flipped_faces;

    // modifiable properties
    bool is_tangent = false;
    bool is_used    = false;
    int base_label  = -1;
    int chain       = -1;
    int curve       = -1;

    // tangent constraints for fitting
    bool has_t0_constraint = false;
    double t0x, t0y;

    bool has_t1_constraint = false;
    double t1x, t1y;

    // adjacency
    std::vector<int> qedges; // list of adjacent qedges
  };

} // namespace IGSV
