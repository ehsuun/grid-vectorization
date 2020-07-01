//=============================================================================
//
//  FUNCTION : get_narrow_band_mesh
//
//=============================================================================

#pragma once

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================

  bool get_narrow_band_mesh(const Eigen::MatrixXd& _V,                                 //
                            const Eigen::MatrixXi& _F,                                 //
                            const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                            Eigen::MatrixXd& _V_narrow,                                //
                            Eigen::MatrixXi& _F_narrow                                 //
  );

} // namespace IGSV
