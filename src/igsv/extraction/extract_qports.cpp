#include <igsv/extraction/extract_qports.h>
#include <igsv/extraction/exact_predicates.h> // typedefs: complexd, complexi, EP_IC_K, Point_2, Vector_2
#include <igsv/extraction/get_dir_as_complex.h>

#include <igsv/extraction_entities/QPort.h>
#include <igsv/extraction_entities/QVertex.h>

#include <igsv/common/defs.h>

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

namespace IGSV {

  using TriCoords = CGAL::Barycentric_coordinates::Triangle_coordinates_2<EP_IC_K>;

  // ===========================================================================

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
                      std::vector<std::pair<int, int>>& _QP_QP_connections) {

    _F_qports.clear();
    _F_qports.resize(_F.rows());

    _QP_QP_connections.clear();

    int qpi = 0;

    int n_V_qports = 0;
    int n_E_qports = 0;
    int n_F_qports = 0;

    // loop over q-vertices
    for (int qvi = 0; qvi < _QV.size(); ++qvi) {

      switch (_QV[qvi].type) {

      // ===========================================================================
      // == Q-PORTS : VERTEX Q-VERTICES ============================================
      // ===========================================================================
      case QVertexType::VERTEX: {

        // find vertex id
        int vid;
        if (_QV[qvi].bc(0) > 0.5)
          vid = _F(_QV[qvi].f, 0);
        else if (_QV[qvi].bc(1) > 0.5)
          vid = _F(_QV[qvi].f, 1);
        else
          vid = _F(_QV[qvi].f, 2);

        // loop over outgoing edges
        for (int i = 0; i < _VT[vid].size(); ++i) {

          const int f = _VT[vid][i];  // face index
          const int k = _VTi[vid][i]; // index of vid in the face

          if (!_F_isNarrow(f))
            continue; // skip if needed

          // xy coords of the triangle
          const complexd x(_V(_F(f, (k + 0) % 3), 0), _V(_F(f, (k + 0) % 3), 1));
          const complexd y(_V(_F(f, (k + 1) % 3), 0), _V(_F(f, (k + 1) % 3), 1));
          const complexd z(_V(_F(f, (k + 2) % 3), 0), _V(_F(f, (k + 2) % 3), 1));

          // uv coords of the triangle
          const complexd u(_UV(_FUV(f, (k + 0) % 3), 0), _UV(_FUV(f, (k + 0) % 3), 1));
          const complexd v(_UV(_FUV(f, (k + 1) % 3), 0), _UV(_FUV(f, (k + 1) % 3), 1));
          const complexd w(_UV(_FUV(f, (k + 2) % 3), 0), _UV(_FUV(f, (k + 2) % 3), 1));

          // barycentric coords computation
          TriCoords tri_uv_coords(Point_2(u.real(), u.imag()), //
                                  Point_2(v.real(), v.imag()), //
                                  Point_2(w.real(), w.imag()));

          // orientation
          const int orient = ORIENT2D(u, v, w);

          // skip if degenerate
          if (orient == 0)
            continue;

          int r = 0;
          if (orient > 0) {
            while (POINTS_INTO(get_dir_as_complex(r), u, v, w))
              r++; // rotate ccw until d points OUTSIDE of f
          } else {
            while (POINTS_INTO(get_dir_as_complex(r), u, w, v))
              r++; // rotate ccw until d points OUTSIDE of f
          }

          for (int j = 1; j < 4; j++) {

            if (orient > 0) { // regular, rotate in the cw direction (negative)
              r--;
              if (r < 0)
                r += 4;
            } else { // flipped, rotate in the ccw direction (positive)
              r++;
            }
            r %= 4;

            const complexi d = get_dir_as_complex(r);

            const bool add_new_port = (orient > 0) ? (POINTS_INTO(d, u, v, w) || IS_COLLINEAR(v - u, d))
                                                   : (POINTS_INTO(d, u, w, v) || IS_COLLINEAR(v - u, d));

            if (add_new_port) {

              // endpoint in uv space
              const complexd u_end = u + 0.5 * complexd(d.real(), d.imag());
              // barycentric coords w.r.t. f
              std::vector<double> bc;
              tri_uv_coords(Point_2(u_end.real(), u_end.imag()), std::inserter(bc, bc.end()));
              const complexd q_end = bc[0] * x + bc[1] * y + bc[2] * z;

              _QP.emplace_back(qvi, r, f, std::round(u.real()), std::round(u.imag()), q_end.real(), q_end.imag());

              _QV[qvi].qports.push_back(qpi);

              _F_qports[f].push_back(qpi);

              qpi++;
              n_V_qports++;
            }
          }
        }
        break;
      }

      // ===========================================================================
      // == Q-PORTS : EDGE Q-VERTICES ==============================================
      // ===========================================================================
      case QVertexType::EDGE: {
        break;
      }

      // ===========================================================================
      // == Q-PORTS : FACE Q-VERTICES ==============================================
      // ===========================================================================
      case QVertexType::FACE: {

        const int f = _QV[qvi].f;

        // xy coords of the triangle
        const complexd x(_V(_F(f, 0), 0), _V(_F(f, 0), 1));
        const complexd y(_V(_F(f, 1), 0), _V(_F(f, 1), 1));
        const complexd z(_V(_F(f, 2), 0), _V(_F(f, 2), 1));

        // uv coords of the triangle
        const complexd u(_UV(_FUV(f, 0), 0), _UV(_FUV(f, 0), 1));
        const complexd v(_UV(_FUV(f, 1), 0), _UV(_FUV(f, 1), 1));
        const complexd w(_UV(_FUV(f, 2), 0), _UV(_FUV(f, 2), 1));

        // barycentric coords computation
        TriCoords tri_uv_coords(Point_2(u.real(), u.imag()), //
                                Point_2(v.real(), v.imag()), //
                                Point_2(w.real(), w.imag()));

        // orientation
        const int orient = ORIENT2D(u, v, w);

        // skip if degenerate
        if (orient == 0)
          continue;

        // uv coords of the q vertex w.r.t. f
        const complexd u_star = _QV[qvi].bc(0) * u + _QV[qvi].bc(1) * v + _QV[qvi].bc(2) * w;
        const complexi u_star_i(std::round(u_star.real()),
                                std::round(u_star.imag())); // IMPORTANT: keep the round!

        // add one port for each direction
        for (int r = 0; r < 4; r++) { // rotate int the cw direction (negative)

          int r_star = r; // flip if orientation is flipped
          if (orient < 0)
            r_star = (r_star + 2) % 4;

          // get the endpoint in xy space, for visu
          const complexi d     = get_dir_as_complex(r_star);
          const complexd u_end = u_star + 0.5 * complexd(d.real(), d.imag());

          std::vector<double> bc;
          tri_uv_coords(Point_2(u_end.real(), u_end.imag()), std::inserter(bc, bc.end()));
          const complexd q_end = bc[0] * x + bc[1] * y + bc[2] * z; // xy coords of u_end

          _QP.emplace_back(qvi, r_star, f, u_star_i.real(), u_star_i.imag(), q_end.real(), q_end.imag());

          _QV[qvi].qports.push_back(qpi);

          _F_qports[f].push_back(qpi);

          qpi++;
          n_F_qports++;
        }

        // add connections between last four
        _QP_QP_connections.emplace_back(qpi - 4, qpi - 2);
        _QP_QP_connections.emplace_back(qpi - 3, qpi - 1);
        break;
      }

      default:
        DEBUG_PRINT_WARNING("unexpected QVertexType in extract_qports");
        return false;
      }

    } // end loop over q-vertices

    DEBUG_PRINT_INFO("%d vertex q-ports", n_V_qports);
    DEBUG_PRINT_INFO("%d edge q-ports", n_E_qports);
    DEBUG_PRINT_INFO("%d face q-ports", n_F_qports);
    DEBUG_PRINT_INFO("%zd total q-ports", _QP.size());
    DEBUG_PRINT_INFO("%zd q-port connections", _QP_QP_connections.size());

    return true;
  }

} // namespace IGSV
