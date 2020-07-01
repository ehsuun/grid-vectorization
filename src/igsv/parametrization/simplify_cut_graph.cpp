#include <igsv/parametrization/simplify_cut_graph.h>
#include <Eigen/Core>
#include <vector>

namespace IGSV {

  // ===========================================================================

  bool simplify_cut_graph_iterate(const Eigen::MatrixXd& _V,                                         //
                                  const Eigen::MatrixXi& _F,                                         //
                                  const Eigen::MatrixXi& _TT,                                        //
                                  const Eigen::MatrixXi& _TTi,                                       //
                                  const Eigen::VectorXi& _V_singIndex,                               //
                                  const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrowBoundary, //
                                  Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_customCut) {
    int n_iter          = 0;
    int n_removed_total = 0;
    int n_removed_iter  = 0;
    do {
      if (n_iter > 1000) { // safety check
        // std::cerr << "WARNING in " << __FUNCTION__ << " : max # of iterations reached\n";
        return false;
      }
      n_removed_iter = simplify_cut_graph(_V, _F, _TT, _TTi, _V_singIndex, _V_isNarrowBoundary, _E_customCut);
      n_removed_total += n_removed_iter;
      n_iter++;
    } while (n_removed_iter > 0);
    // std::cout << "  removed " << n_removed_total << " vertices in " << n_iter << " iterations\n";
    return true;
  }

  // ===========================================================================

  int simplify_cut_graph(const Eigen::MatrixXd& _V,                                         //
                         const Eigen::MatrixXi& _F,                                         //
                         const Eigen::MatrixXi& _TT,                                        //
                         const Eigen::MatrixXi& _TTi,                                       //
                         const Eigen::VectorXi& _V_singIndex,                               //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrowBoundary, //
                         Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_customCut) {

    const int nv = _V.rows();
    const int nf = _F.rows();

    int n_removed_iter = 0;

    // the edges that are in the cut graph are stored in _E_customCut
    // find and remove edges with valence 1 vertices (that are not on the boundary)

    // what is the valence of vertices in the current cut graph?
    std::vector<int> V_valence(nv, 0);

    // list of cut edges for each vertex
    std::vector<std::vector<std::pair<int, int>>> V_cut_edges(nv);

    Eigen::Matrix<bool, Eigen::Dynamic, 3> e_visited;
    e_visited.resize(nf, 3);
    e_visited.fill(false);

    for (int f0 = 0; f0 < nf; ++f0)
      for (int k0 = 0; k0 < 3; ++k0)

        if (!e_visited(f0, k0) && _E_customCut(f0, k0)) {

          const int v0 = _F(f0, k0);
          const int v1 = _F(f0, (k0 + 1) % 3);

          V_valence[v0]++;
          V_valence[v1]++;

          V_cut_edges[v0].push_back(std::make_pair(f0, k0));
          V_cut_edges[v1].push_back(std::make_pair(f0, k0));

          e_visited(f0, k0) = true;

          const int f1 = _TT(f0, k0);
          const int k1 = _TTi(f0, k0);

          if (f1 != -1)
            e_visited(f1, k1) = true;
        }

    // find valence-1 vertices inside the narrow band
    for (int v = 0; v < nv; ++v)
      if (!_V_isNarrowBoundary[v] && (V_valence[v] == 1) && (_V_singIndex(v) == 0)) {

        // remove the edges
        for (const auto& edge : V_cut_edges[v]) {
          const int f = edge.first;
          const int k = edge.second;

          _E_customCut(f, k) = false;
          if (_TT(f, k) != -1)
            _E_customCut(_TT(f, k), _TTi(f, k)) = false;
        }

        // remove the vertex
        n_removed_iter++;
      }

    return n_removed_iter;
  }

  // ===========================================================================

} // namespace IGSV
