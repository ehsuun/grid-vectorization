//=============================================================================
//
//  FUNCTION : extract_qports
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
  struct QPort;

  //== FUNCTION DEFINITION ====================================================

  bool extract_qports(const Eigen::MatrixXd& _V,                                 //
                      const Eigen::MatrixXi& _F,                                 //
                      const Eigen::MatrixXd& _UV,                                //
                      const Eigen::MatrixXi& _FUV,                               //
                      const std::vector<std::vector<int>>& _VT,                  //
                      const std::vector<std::vector<int>>& _VTi,                 //
                      const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                      std::vector<QVertex>& _QV,                                 //
                      std::vector<QPort>& _QP,                                   //
                      std::vector<std::vector<int>>& _F_qports,                  //
                      std::vector<std::pair<int, int>>& _QP_QP_connections);

} // namespace IGSV
