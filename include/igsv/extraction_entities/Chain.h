//=============================================================================
//
//  STRUCT : Chain
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <map>
#include <set>
#include <vector>

#include <Eigen/Core>

#include "Sample.h"

//== NAMESPACES ===============================================================

namespace IGSV {

  //== STRUCT DEFINITION ======================================================

  struct Chain {
    int qvi0; // source qvertex
    int qvi1; // target qvertex

    // std::map<std::pair<int, int>, int> qverts2qedge;

    std::set<int> qedges;          // list of qedges
    std::set<int> adjacent_qedges; // list of qedges adjacent to this chain

    std::vector<int> qverts;        // list of qvertices
    std::vector<int> qverts_degree; // degrees of qvertices
    std::map<int, int> vindex_map;  // map from global index to local index
    std::vector<double> t;          // harmonic parametrization of qvertices
    std::vector<ChainSample> samples;

    int curve = -1;

    Eigen::VectorXi histogram;
    bool is_removed = false;

    /*
      http://www.alglib.net/translator/man/manual.cpp.html#sub_spline1dunpack

      coefficients table, unpacked format, array[0..N-2, 0..5].
        For I = 0...N-2:
            Tbl[I,0] = X[i]
            Tbl[I,1] = X[i+1]
            Tbl[I,2] = C0
            Tbl[I,3] = C1
            Tbl[I,4] = C2
            Tbl[I,5] = C3
        On [x[i], x[i+1]] spline is equals to:
            S(x) = C0 + C1*t + C2*t^2 + C3*t^3
            t = x-x[i]
    */
    Eigen::MatrixXd CubicSpline;

    // 4 x n_segments matrix, each column: control polygon of a cubic Bezier curve
    Eigen::MatrixXcd CP_spline;

    double fitting_cost    = 1.e6;
    double normalized_cost = 1.e6;
  };

} // namespace IGSV
