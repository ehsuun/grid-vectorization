//=============================================================================
//
//  FUNCTION : extract_qvertices
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  struct QVertex;

  //== FUNCTION DEFINITION ====================================================

  bool extract_qvertices(const Eigen::MatrixXd& _V,                                 //
                         const Eigen::MatrixXi& _F,                                 //
                         const Eigen::MatrixXd& _UV,                                //
                         const Eigen::MatrixXi& _FUV,                               //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                         std::vector<QVertex>& _QV);

} // namespace IGSV
