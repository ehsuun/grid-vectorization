#pragma once

#include <Eigen/Core>

namespace IGSV {

  int simplify_cut_graph(                                                //
      const Eigen::MatrixXd& _V,                                         //
      const Eigen::MatrixXi& _F,                                         //
      const Eigen::MatrixXi& _TT,                                        //
      const Eigen::MatrixXi& _TTi,                                       //
      const Eigen::VectorXi& _V_singIndex,                               //
      const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrowBoundary, //
      Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_customCut);

  bool simplify_cut_graph_iterate(                                       //
      const Eigen::MatrixXd& _V,                                         //
      const Eigen::MatrixXi& _F,                                         //
      const Eigen::MatrixXi& _TT,                                        //
      const Eigen::MatrixXi& _TTi,                                       //
      const Eigen::VectorXi& _V_singIndex,                               //
      const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrowBoundary, //
      Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_customCut);

} // namespace IGSV
