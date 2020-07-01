//=============================================================================
//
//  FUNCTION : find_nearest
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include "_cgal_delaunay_2.h" // using IndexedDT
#include <Eigen/Core>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================

  bool find_nearest(const Eigen::MatrixXd& _V,                //
                    const Eigen::MatrixXi& _F,                //
                    const std::vector<std::vector<int>>& _VT, //
                    const Eigen::VectorXcd& _PX_XY,           //
                    IndexedDT& _tri,                          //
                    Eigen::VectorXi& _PX_nearestID,           //
                    Eigen::VectorXcd& _PX_nearestXY,          //
                    Eigen::VectorXi& _PX_fid,                 //
                    Eigen::MatrixXd& _PX_bc);

} // namespace IGSV
