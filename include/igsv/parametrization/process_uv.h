//=============================================================================
//
//  FUNCTION :  sanitize_parametrization
//              precompute_transition_fn
//
//=============================================================================

#pragma once

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================

  bool sanitize_parametrization(Eigen::MatrixXd& _UV);

  bool precompute_transition_fn(const Eigen::MatrixXd& _UV,            //
                                const Eigen::MatrixXi& _FUV,           //
                                const Eigen::MatrixXi& _TT,            //
                                const Eigen::MatrixXi& _TTi,           //
                                const Eigen::MatrixXi& _E_periodJumps, //
                                Eigen::MatrixXcd& _R_transfn,          //
                                Eigen::MatrixXcd& _T_transfn,          //
                                Eigen::MatrixXcd& _T_integer_transfn,  //
                                Eigen::MatrixXcd& _T_decimal_transfn,  //
                                Eigen::Matrix<bool, Eigen::Dynamic, 3>& _flip_transfn);

} // namespace IGSV
