#include <igsv/extraction/parametrize_chains.h>
#include <igsv/common/defs.h>
#include <igsv/extraction_entities/Chain.h>
#include <igsv/extraction_entities/QEdge.h>
#include <igsv/extraction_entities/QVertex.h>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <igl/slice.h>

namespace IGSV {

  using TripletD = Eigen::Triplet<double>;

  //===========================================================================

  bool parametrize_chains(const std::vector<QVertex>& _QV, //
                          const std::vector<QEdge>& _QE,   //
                          std::vector<Chain>& _CH) {

    int chid = 0;
    for (auto& _chain : _CH) {

      if (_chain.qvi0 == _chain.qvi1) {

        const int n_qverts = _chain.qverts.size() + 1;
        const int n_qedges = _chain.qedges.size();

        if (n_qedges == 0)
          return false;

        //// geometric Laplacian (inverse lengths as edge weights)
        std::vector<TripletD> ijw;
        ijw.reserve(n_qedges * 4);
        int i0, i1;
        bool first = true;
        double w;
        for (auto qei : _chain.qedges) {

          if (_QE[qei].qvi0 == _chain.qvi0) {
            if (first) {
              i0    = 0;
              first = false;
            } else {
              i0 = n_qverts - 1;
            }
          } else {
            i0 = _chain.vindex_map[_QE[qei].qvi0];
          }

          if (_QE[qei].qvi1 == _chain.qvi0) {
            if (first) {
              i1    = 0;
              first = false;
            } else {
              i1 = n_qverts - 1;
            }
          } else {
            i1 = _chain.vindex_map[_QE[qei].qvi1];
          }

          w = 1. / _QE[qei].length;
          ijw.push_back(TripletD(i0, i0, +w));
          ijw.push_back(TripletD(i1, i0, -w));
          ijw.push_back(TripletD(i1, i1, +w));
          ijw.push_back(TripletD(i0, i1, -w));
        }
        Eigen::SparseMatrix<double> L;
        L.resize(n_qverts, n_qverts);
        L.setFromTriplets(ijw.begin(), ijw.end());

        //// fixed vertices : first and last
        const int i_source    = 0;
        const int i_target    = n_qverts - 1;
        const double t_source = 0.0;
        const double t_target = 1.0;

        //// indices of fixed
        Eigen::VectorXi idx_fixd(2);
        idx_fixd << i_source, i_target;

        //// indices of free
        Eigen::VectorXi idx_free(n_qverts - 2);
        int row = 0;
        for (int i = 0; i < n_qverts; ++i) {
          if (i == i_source)
            continue;
          if (i == i_target)
            continue;
          idx_free(row++) = i;
        }

        //// slice matrices
        Eigen::SparseMatrix<double> Lc, Lf;
        igl::slice(L, idx_free, idx_free, Lf);
        igl::slice(L, idx_free, idx_fixd, Lc);

        //// rhs
        Eigen::VectorXd Bc(2);
        Bc << t_source, t_target;

        //// solve via Cholesky
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_param(Lf);
        if (solver_param.info() != Eigen::Success) {
          DEBUG_PRINT_WARNING("chain #%d: fail in Eigen::LDLT decomposition", chid);
          return false;
        }

        const Eigen::VectorXd Xf = solver_param.solve(-Lc * Bc);
        if (solver_param.info() != Eigen::Success) {
          DEBUG_PRINT_WARNING("chain #%d: fail in Eigen::LDLT solve", chid);
          return false;
        }

        //// copy the result
        Eigen::VectorXd X(n_qverts);
        X(i_source) = t_source;
        X(i_target) = t_target;
        for (int i = 0; i < n_qverts - 2; ++i)
          X(idx_free(i)) = Xf(i);

        _chain.t.resize(X.size());
        Eigen::VectorXd::Map(&(_chain.t[0]), X.size()) = X;

      } else {

        const int n_qverts = _chain.qverts.size();
        const int n_qedges = _chain.qedges.size();

        if (n_qedges == 0)
          return false;

        //// geometric Laplacian (inverse lengths as edge weights)
        std::vector<TripletD> ijw;
        ijw.reserve(n_qedges * 4);
        int i0, i1;
        double w;
        for (auto qei : _chain.qedges) {
          i0 = _chain.vindex_map[_QE[qei].qvi0];
          i1 = _chain.vindex_map[_QE[qei].qvi1];
          w  = 1. / _QE[qei].length;
          ijw.push_back(TripletD(i0, i0, +w));
          ijw.push_back(TripletD(i1, i0, -w));
          ijw.push_back(TripletD(i1, i1, +w));
          ijw.push_back(TripletD(i0, i1, -w));
        }
        Eigen::SparseMatrix<double> L;
        L.resize(n_qverts, n_qverts);
        L.setFromTriplets(ijw.begin(), ijw.end());

        //// fixed vertices : first and last
        const int i_source    = _chain.vindex_map[_chain.qvi0];
        const int i_target    = _chain.vindex_map[_chain.qvi1];
        const double t_source = 0.0;
        const double t_target = 1.0;

        //// indices of fixed
        Eigen::VectorXi idx_fixd(2);
        idx_fixd << i_source, i_target;

        //// indices of free
        Eigen::VectorXi idx_free(n_qverts - 2);
        int row = 0;
        for (int i = 0; i < n_qverts; ++i) {
          if (i == i_source)
            continue;
          if (i == i_target)
            continue;
          idx_free(row++) = i;
        }

        //// slice matrices
        Eigen::SparseMatrix<double> Lc, Lf;
        igl::slice(L, idx_free, idx_free, Lf);
        igl::slice(L, idx_free, idx_fixd, Lc);

        //// rhs
        Eigen::VectorXd Bc(2);
        Bc << t_source, t_target;

        //// solve via Cholesky
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver_param(Lf);
        if (solver_param.info() != Eigen::Success) {
          DEBUG_PRINT_WARNING("chain #%d: fail in Eigen::LDLT decomposition", chid);
          return false;
        }

        const Eigen::VectorXd Xf = solver_param.solve(-Lc * Bc);
        if (solver_param.info() != Eigen::Success) {
          DEBUG_PRINT_WARNING("chain #%d: fail in Eigen::LDLT solve", chid);
          return false;
        }

        //// copy the result
        Eigen::VectorXd X(n_qverts);
        X(i_source) = t_source;
        X(i_target) = t_target;
        for (int i = 0; i < n_qverts - 2; ++i)
          X(idx_free(i)) = Xf(i);

        _chain.t.resize(X.size());
        Eigen::VectorXd::Map(&(_chain.t[0]), X.size()) = X;
      }

      chid++;
    }

    return true;
  }

} // namespace IGSV
