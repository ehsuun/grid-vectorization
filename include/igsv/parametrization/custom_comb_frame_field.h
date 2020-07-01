//=============================================================================
//
//  FUNCTION : custom_comb_frame_field
//
//=============================================================================

#pragma once

//== NAMESPACES ===============================================================

namespace IGSV {

  //== FUNCTION DEFINITION ====================================================

  bool custom_comb_frame_field(const Eigen::MatrixXd& _V,                                  //
                               const Eigen::MatrixXi& _F,                                  //
                               const Eigen::MatrixXi& _TT,                                 //
                               const Eigen::MatrixXi& _TTi,                                //
                               const std::vector<std::vector<int>>& _VT,                   //
                               const std::vector<std::vector<int>>& _VTi,                  //
                               const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrow,  //
                               const Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrow,  //
                               const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow,  //
                               const Eigen::VectorXd& _X0_arg,                             //
                               const Eigen::VectorXd& _X1_arg,                             //
                               const Eigen::VectorXd& _X0_mag,                             //
                               const Eigen::VectorXd& _X1_mag,                             //
                               Eigen::VectorXi& _V_singIndex,                              //
                               Eigen::MatrixXi& _E_periodJumps,                            //
                               Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_customCut,       //
                               std::vector<std::pair<int, int>>& _TT_dualConnections,      //
                               Eigen::MatrixXd& _BIS0_unit,                                //
                               Eigen::MatrixXd& _BIS1_unit,                                //
                               Eigen::MatrixXd& _BIS0_unit_comb,                           //
                               Eigen::MatrixXd& _BIS1_unit_comb,                           //
                               Eigen::MatrixXd& _X0_unit,                                  //
                               Eigen::MatrixXd& _X1_unit,                                  //
                               Eigen::MatrixXd& _X0_nonunit,                               //
                               Eigen::MatrixXd& _X1_nonunit,                               //
                               Eigen::MatrixXd& _X0_unit_comb,                             //
                               Eigen::MatrixXd& _X1_unit_comb,                             //
                               Eigen::MatrixXd& _X0_nonunit_comb,                          //
                               Eigen::MatrixXd& _X1_nonunit_comb,                          //
                               Eigen::VectorXcd (&_F_frameField)[4]);

} // namespace IGSV
