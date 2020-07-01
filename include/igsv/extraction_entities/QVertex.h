//=============================================================================
//
//  STRUCT : QVertex
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== ENUMS ==================================================================
  enum QVertexType { VERTEX = 0, EDGE, FACE };

  //== STRUCT DEFINITION ======================================================
  struct QVertex {

    QVertex(double _x, double _y, int _f, const Eigen::Vector3d& _bc, QVertexType _type)
        : x(_x), y(_y), f(_f), bc(_bc), type(_type) {}

    // const properties
    const double x;           // x-coord in the sketch space
    const double y;           // y-coord in the sketch space
    const int f;              // index of face where the qvertex is located (one of the faces if type==EDGE or FACE)
    const Eigen::Vector3d bc; // barycentric coords w.r.t. f
    const QVertexType type;

    // modifiable properties
    bool is_split_node = false; // should chains be split at this qvertex?
    int valence        = 0;     // what is the valence of this qvertex? (w.r.t. active qedges)

    // adjacency
    std::vector<int> qports; // list of associated qports
    std::vector<int> qedges; // list of associated qedges
  };

} // namespace IGSV
