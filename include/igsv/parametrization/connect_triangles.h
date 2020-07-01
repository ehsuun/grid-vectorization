//=============================================================================
//
//  FUNCTION : connect_triangles
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <opencv2/opencv.hpp>
#include <vector>

#include "_cgal_delaunay_2.h" // using IndexedDT

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================

  void connect_triangles(const Eigen::MatrixXd& _V,  //
                         const Eigen::MatrixXi& _F,  //
                         const Eigen::MatrixXi& _TT, //
                         const IndexedDT& _tri,      //
                         const int f0,               //
                         const int f1,               //
                         std::vector<int>& faceList, //
                         std::vector<int>& edgeList);

  void connect_triangles(const Eigen::MatrixXi& _TT, //
                         const IndexedDT& _tri,      //
                         const int f0,               //
                         const int f1,               //
                         std::vector<int>& faceList, //
                         std::vector<int>& edgeList, //
                         const double p0_x,          //
                         const double p0_y,          //
                         const double p1_x,          //
                         const double p1_y,          //
                         const IndexedDT::Face_handle fh0 = IndexedDT::Face_handle());

}; // namespace IGSV
