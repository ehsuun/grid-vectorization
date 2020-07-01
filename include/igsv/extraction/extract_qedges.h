//=============================================================================
//
//  FUNCTION : extract_qedges
//
//=============================================================================

#pragma once

//== INCLUDES =================================================================

#include <Eigen/Core>
#include <map>
#include <vector>

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FORWARD DECLARATIONS ===================================================

  struct QVertex;
  struct QPort;
  struct QEdge;

  //== FUNCTION DEFINITION ====================================================

  bool extract_qedges(const Eigen::MatrixXd& _UV,                                //
                      const Eigen::MatrixXi& _FUV,                               //
                      const Eigen::MatrixXi& _TT,                                //
                      const Eigen::MatrixXi& _TTi,                               //
                      const Eigen::MatrixXi& _E_periodJumps,                     //
                      const Eigen::MatrixXcd& _R_transfn,                        //
                      const Eigen::MatrixXcd& _T_integer_transfn,                //
                      const Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrow, //
                      const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                      const Eigen::VectorXi& _F_labels,                          //
                      const std::vector<std::vector<int>>& _F_qports,            //
                      std::vector<QVertex>& _QV,                                 //
                      std::vector<QPort>& _QP,                                   //
                      std::vector<QEdge>& _QE,                                   //
                      std::vector<std::pair<int, int>>& _QP_QP_connections,      //
                      std::map<std::pair<int, int>, int>& _QVPair_to_QE          //
  );

} // namespace IGSV
