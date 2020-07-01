//=============================================================================
//
//  FUNCTION : store_orientation
//
//=============================================================================

#pragma once

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================

  bool store_orientation(const Eigen::MatrixXd& _UV,                                //
                         const Eigen::MatrixXi& _FUV,                               //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                         Eigen::VectorXi& _F_uvOrientation);

} // namespace IGSV
