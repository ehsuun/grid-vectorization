#include <igsv/extraction/extract_qvertices.h>

#include <igsv/common/defs.h>
#include <igsv/extraction_entities/QVertex.h>

#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>
#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_2.h>

using EP_IC_K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using TriCoords = CGAL::Barycentric_coordinates::Triangle_coordinates_2<EP_IC_K>;
using SegCoords = CGAL::Barycentric_coordinates::Segment_coordinates_2<EP_IC_K>;
using Point_2   = EP_IC_K::Point_2;

namespace IGSV {

  // ===========================================================================

  bool extract_VERTEX_qvertices(const Eigen::MatrixXd& _V,                                 //
                                const Eigen::MatrixXi& _F,                                 //
                                const Eigen::MatrixXd& _UV,                                //
                                const Eigen::MatrixXi& _FUV,                               //
                                const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                                std::vector<QVertex>& _QV) {

    // ===========================================================================
    // == VERTEX Q-VERTICES ======================================================
    // ===========================================================================

    int n_V_qverts = 0;
    Eigen::Vector3d bc;

    Eigen::Matrix<bool, Eigen::Dynamic, 1> v_done(_V.rows());
    v_done.fill(false);

    for (int f = 0; f < _F.rows(); ++f) { // loop over faces

      if (!_F_isNarrow(f))
        continue; // skip white faces

      for (int k = 0; k < 3; ++k) { // loop over corners

        const int ci = _FUV(f, k); // get indices
        const int vi = _F(f, k);   //

        if (v_done(vi))
          continue; // skip if vertex already processed

        v_done(vi) = true; // mark vertex as processed

        const double u = _UV(ci, 0); // uv coords of the corner
        const double v = _UV(ci, 1); //

        if (std::abs(u - std::round(u)) < 1e-8 && // check if u and v are integers
            std::abs(v - std::round(v)) < 1e-8) {
          bc.fill(0);
          bc(k) = 1.;
          _QV.emplace_back(_V(vi, 0), _V(vi, 1), f, bc, QVertexType::VERTEX);
          n_V_qverts++;
        }
      }
    }

    DEBUG_PRINT_INFO("%d vertex q-vertices", n_V_qverts);

    return true;
  }

  // ===========================================================================

  bool extract_EDGE_qvertices(const Eigen::MatrixXd& _V,                                 //
                              const Eigen::MatrixXi& _F,                                 //
                              const Eigen::MatrixXd& _UV,                                //
                              const Eigen::MatrixXi& _FUV,                               //
                              const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                              std::vector<QVertex>& _QV) {

    // ===========================================================================
    // == EDGE Q-VERTICES ========================================================
    // ===========================================================================

    int n_E_qverts = 0;
    Eigen::Vector3d bc;
    Eigen::Matrix<double, 3, 2> Tuv;
    double u_lower, u_upper, v_lower, v_upper;

    for (int f = 0; f < _F.rows(); ++f) { // loop over faces

      if (!_F_isNarrow(f))
        continue; // skip white faces

      // get integer bounds
      for (int k = 0; k < 3; k++)
        Tuv.row(k) << _UV.row(_FUV(f, k));

      u_lower = std::ceil(Tuv.col(0).minCoeff());
      u_upper = std::floor(Tuv.col(0).maxCoeff()) + 1;
      v_lower = std::ceil(Tuv.col(1).minCoeff());
      v_upper = std::floor(Tuv.col(1).maxCoeff()) + 1;

      for (int k = 0; k < 3; ++k) { // loop over edges

        const int v0 = _F(f, k);
        const int v1 = _F(f, (k + 1) % 3);

        if (v0 > v1)
          continue; // only process one oriented edge, the one where v0 < v1

        const int c0 = _FUV(f, k);
        const int c1 = _FUV(f, (k + 1) % 3);

        const EP_IC_K::Segment_2 seg_uv(Point_2(_UV(c0, 0), _UV(c0, 1)), //
                                        Point_2(_UV(c1, 0), _UV(c1, 1)));

        const EP_IC_K::Segment_2 seg_xy(Point_2(_V(v0, 0), _V(v0, 1)), //
                                        Point_2(_V(v1, 0), _V(v1, 1)));

        SegCoords seg_uv_coords(seg_uv[0], seg_uv[1]);

        // TODO : too many tests, modify
        for (int u = u_lower; u < u_upper; ++u) {
          for (int v = v_lower; v < v_upper; ++v) {

            double pu = u;
            double pv = v;

            Point_2 p_uv(pu, pv);

            if (seg_uv.has_on(p_uv) && (seg_uv[0] != p_uv) && (seg_uv[1] != p_uv)) {

              // the integer point is inside the edge segment

              // get barycentric coords
              std::vector<double> bc_vec;
              seg_uv_coords(p_uv, std::inserter(bc_vec, bc_vec.end()));

              // compute the position in xy space
              double px = 0, py = 0;
              for (int k = 0; k < 2; ++k) {
                px += bc_vec[k] * seg_xy[k].x();
                py += bc_vec[k] * seg_xy[k].y();
              }

              bc(k)           = bc_vec[0];
              bc((k + 1) % 3) = bc_vec[1];
              bc((k + 2) % 3) = 0;
              _QV.emplace_back(px, py, f, bc, QVertexType::EDGE);

              n_E_qverts++;

            } // if point on segment
          }   // loop over v
        }     // loop over u
      }       // loop over edges
    }         // loop over faces

    DEBUG_PRINT_INFO("%d edge q-vertices", n_E_qverts);

    return true;
  }

  // ===========================================================================

  bool extract_FACE_qvertices(const Eigen::MatrixXd& _V,                                 //
                              const Eigen::MatrixXi& _F,                                 //
                              const Eigen::MatrixXd& _UV,                                //
                              const Eigen::MatrixXi& _FUV,                               //
                              const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                              std::vector<QVertex>& _QV) {

    // ===========================================================================
    // == FACE Q-VERTICES ========================================================
    // ===========================================================================

    int n_F_qverts = 0;
    Eigen::Vector3d bc;
    Eigen::Matrix<double, 3, 2> Tuv;
    double u_lower, u_upper, v_lower, v_upper;

    for (int f = 0; f < _F.rows(); ++f) {

      // skip white faces
      if (!_F_isNarrow(f))
        continue;

      // get integer bounds for this triangle
      for (int k = 0; k < 3; k++)
        Tuv.row(k) << _UV.row(_FUV(f, k));

      u_lower = std::ceil(Tuv.col(0).minCoeff());
      u_upper = std::floor(Tuv.col(0).maxCoeff()) + 1;
      v_lower = std::ceil(Tuv.col(1).minCoeff());
      v_upper = std::floor(Tuv.col(1).maxCoeff()) + 1;

      // CGAL : triangles
      // xy space
      const int v0 = _F(f, 0);
      const int v1 = _F(f, 1);
      const int v2 = _F(f, 2);
      const EP_IC_K::Triangle_2 tri_xy(Point_2(_V(v0, 0), _V(v0, 1)), //
                                       Point_2(_V(v1, 0), _V(v1, 1)), //
                                       Point_2(_V(v2, 0), _V(v2, 1)));

      // uv space
      const int c0 = _FUV(f, 0);
      const int c1 = _FUV(f, 1);
      const int c2 = _FUV(f, 2);
      const EP_IC_K::Triangle_2 tri_uv(Point_2(_UV(c0, 0), _UV(c0, 1)), //
                                       Point_2(_UV(c1, 0), _UV(c1, 1)), //
                                       Point_2(_UV(c2, 0), _UV(c2, 1)));

      TriCoords tri_uv_coords(tri_uv[0], tri_uv[1], tri_uv[2]);

      for (int u = u_lower; u < u_upper; ++u) {
        for (int v = v_lower; v < v_upper; ++v) {

          double pu = u;
          double pv = v;

          Point_2 p_uv(pu, pv);

          const bool point_in_f = (tri_uv.orientation() == +1 && tri_uv.has_on_positive_side(p_uv)) ||
                                  (tri_uv.orientation() == -1 && tri_uv.has_on_negative_side(p_uv));

          if (point_in_f) {

            // the integer point is inside the triangle

            // get barycentric coords
            std::vector<double> bc_vec;
            tri_uv_coords(p_uv, std::inserter(bc_vec, bc_vec.end()));

            // compute the position in xy space
            double px = 0, py = 0;
            for (int k = 0; k < 3; ++k) {
              px += bc_vec[k] * tri_xy[k].x();
              py += bc_vec[k] * tri_xy[k].y();
            }

            bc(0) = bc_vec[0];
            bc(1) = bc_vec[1];
            bc(2) = bc_vec[2];
            _QV.emplace_back(px, py, f, bc, QVertexType::FACE);

            n_F_qverts++;

          } // test if inside
        }   // loop over v
      }     // loop over u
    }       // loop over faces

    DEBUG_PRINT_INFO("%d face q-vertices", n_F_qverts);

    return true;
  }

  // ===========================================================================

  bool extract_qvertices(const Eigen::MatrixXd& _V,                                 //
                         const Eigen::MatrixXi& _F,                                 //
                         const Eigen::MatrixXd& _UV,                                //
                         const Eigen::MatrixXi& _FUV,                               //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                         std::vector<QVertex>& _QV) {

    if (!extract_VERTEX_qvertices(_V, _F, _UV, _FUV, _F_isNarrow, _QV))
      return false;

    if (!extract_EDGE_qvertices(_V, _F, _UV, _FUV, _F_isNarrow, _QV))
      return false;

    if (!extract_FACE_qvertices(_V, _F, _UV, _FUV, _F_isNarrow, _QV))
      return false;

    DEBUG_PRINT_INFO("%zd total q-vertices", _QV.size());

    return true;
  }

} // namespace IGSV
