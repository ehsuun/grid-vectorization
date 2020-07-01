#include <igsv/extraction/extract_qedges.h>
#include <igsv/extraction/exact_predicates.h> // typedefs: complexd, complexi, EP_IC_K, Point_2, Vector_2
#include <igsv/extraction/get_dir_as_complex.h>

#include <igsv/extraction_entities/QEdge.h>
#include <igsv/extraction_entities/QPort.h>
#include <igsv/extraction_entities/QVertex.h>

#include <igsv/common/CornerType.h>
#include <igsv/common/defs.h>

#include <CGAL/intersections.h>

namespace IGSV {

  using MatrixXb  = Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>;
  using MatrixX3b = Eigen::Matrix<bool, Eigen::Dynamic, 3>;
  using VectorXb  = Eigen::Matrix<bool, Eigen::Dynamic, 1>;

  using Segment_2  = EP_IC_K::Segment_2;
  using Triangle_2 = EP_IC_K::Triangle_2;

  bool extract_qedges(const Eigen::MatrixXd& _UV,                           //
                      const Eigen::MatrixXi& _FUV,                          //
                      const Eigen::MatrixXi& _TT,                           //
                      const Eigen::MatrixXi& _TTi,                          //
                      const Eigen::MatrixXi& _E_periodJumps,                //
                      const Eigen::MatrixXcd& _R_transfn,                   //
                      const Eigen::MatrixXcd& _T_integer_transfn,           //
                      const MatrixX3b& _E_isNarrow,                         //
                      const VectorXb& _F_isNarrow,                          //
                      const Eigen::VectorXi& _F_labels,                     //
                      const std::vector<std::vector<int>>& _F_qports,       //
                      std::vector<QVertex>& _QV,                            //
                      std::vector<QPort>& _QP,                              //
                      std::vector<QEdge>& _QE,                              //
                      std::vector<std::pair<int, int>>& _QP_QP_connections, //
                      std::map<std::pair<int, int>, int>& _QVPair_to_QE) {

    const int n_qverts = _QV.size();
    const int n_qports = _QP.size();

    VectorXb QP_connected;
    QP_connected.resize(n_qports, 1);
    QP_connected.fill(false);

    MatrixXb QV_QV_connections;
    QV_QV_connections.resize(n_qverts, n_qverts);
    QV_QV_connections.fill(false);

    int qei = 0;

    for (int qpi = 0; qpi < n_qports; ++qpi) {

      // const bool SHOW_QEDGE_DEBUG_INFO = qpi == 13 || qpi == 113;
      const bool SHOW_QEDGE_DEBUG_INFO = false;

      if (SHOW_QEDGE_DEBUG_INFO) {
        DEBUG_PRINT_INFO("QPort %d", qpi);
        DEBUG_PRINT_INFO("  f: %d", _QP[qpi].f);
        DEBUG_PRINT_INFO("  r: %d", _QP[qpi].pj);
      }

      if (QP_connected(qpi)) {
        if (SHOW_QEDGE_DEBUG_INFO) {
          DEBUG_PRINT_INFO(" ... skip, already connected.");
        }
        continue;
      }

      const int qvi = _QP[qpi].qvi;
      const int f   = _QP[qpi].f;
      const int r   = _QP[qpi].pj;
      const complexd uv(_QP[qpi].u, _QP[qpi].v);

      if (SHOW_QEDGE_DEBUG_INFO) {
        DEBUG_PRINT_INFO("qpi #%d", qpi);
        DEBUG_PRINT_INFO("    qvi = %d", qvi);
        DEBUG_PRINT_INFO("    f   = %d", f);
        DEBUG_PRINT_INFO("    r   = %d", r);
        DEBUG_PRINT_INFO("    u   = %d", _QP[qpi].u);
        DEBUG_PRINT_INFO("    v   = %d", _QP[qpi].v);
      }

      // init
      int r0 = r;
      int f0 = f;
      int k0 = -1;
      int o_saved;
      complexi a = uv;
      complexi d = get_dir_as_complex(r0);
      complexi b = a + d;

      bool point_is_inside_or_on_boundary = false;
      bool edge_is_tangent                = false;
      bool edge_has_flipped_faces         = false;

      const int maxiter = 1000;
      int iter          = 0;
      while (++iter < maxiter) {

        const int flabel = _F_labels(f0);

        if (SHOW_QEDGE_DEBUG_INFO) {
          DEBUG_PRINT_INFO("  flabel= %d, r0= %d", flabel, r0);
        }

        if ((flabel == CORNER_TYPE_INTEGER_U & r0 % 2 == 0) || (flabel == CORNER_TYPE_INTEGER_V & r0 % 2 == 1)) {
          edge_is_tangent = true;
        }

        // uv coords of the triangle
        const complexd u(_UV(_FUV(f0, 0), 0), _UV(_FUV(f0, 0), 1));
        const complexd v(_UV(_FUV(f0, 1), 0), _UV(_FUV(f0, 1), 1));
        const complexd w(_UV(_FUV(f0, 2), 0), _UV(_FUV(f0, 2), 1));
        const Triangle_2 t(Point_2(u.real(), u.imag()), //
                           Point_2(v.real(), v.imag()), //
                           Point_2(w.real(), w.imag()));

        const int o_current    = t.orientation();
        edge_has_flipped_faces = edge_has_flipped_faces | (o_current < 1);

        if (k0 == -1) {
          o_saved = o_current;
        } else if (o_current != 0 && o_current != o_saved) {
          // change direction and orientation
          o_saved = o_current;
          std::swap(a, b);
          r0 = (r0 + 2) % 4;
          if (SHOW_QEDGE_DEBUG_INFO) {
            DEBUG_PRINT_INFO("  flip : r0 = %d", r0);
          }
        }

        // check if b is inside the triangle f
        const Point_2 pa(a.real(), a.imag());
        const Point_2 pb(b.real(), b.imag());

        point_is_inside_or_on_boundary = t.has_on_bounded_side(pb) | t.has_on_boundary(pb);

        if (SHOW_QEDGE_DEBUG_INFO) {
          DEBUG_PRINT_INFO("  f= %d", f0);
          DEBUG_PRINT_INFO("    u : (%+0.16f, %+0.16f)", t[0].x(), t[0].y());
          DEBUG_PRINT_INFO("    v : (%+0.16f, %+0.16f)", t[1].x(), t[1].y());
          DEBUG_PRINT_INFO("    w : (%+0.16f, %+0.16f)", t[2].x(), t[2].y());
          DEBUG_PRINT_INFO("    a : (%+0.16f, %+0.16f)", pa.x(), pa.y());
          DEBUG_PRINT_INFO("    b : (%+0.16f, %+0.16f)", pb.x(), pb.y());
          DEBUG_PRINT_INFO("    t.has_on_bounded_side(pb)?  %s", t.has_on_bounded_side(pb) ? "YES" : "NO");
          DEBUG_PRINT_INFO("    t.has_on_boundary(pb)?      %s", t.has_on_boundary(pb) ? "YES" : "NO");
          DEBUG_PRINT_INFO("    b == u? %s", pb == t[0] ? "YES" : "NO");
          DEBUG_PRINT_INFO("    b == v? %s", pb == t[1] ? "YES" : "NO");
          DEBUG_PRINT_INFO("    b == w? %s", pb == t[2] ? "YES" : "NO");
        }

        if (point_is_inside_or_on_boundary)
          break;

        // intersect edges with the segment [a,b]
        const Segment_2 ab(pa, pb);
        const Segment_2 e0(t[0], t[1]);
        const Segment_2 e1(t[1], t[2]);
        const Segment_2 e2(t[2], t[0]);
        const auto e0ab = CGAL::intersection(e0, ab);
        const auto e1ab = CGAL::intersection(e1, ab);
        const auto e2ab = CGAL::intersection(e2, ab);

        // get the next edge
        int k0_next = -1;
        if (!e0.has_on(pa) & (k0 != 0))
          if (e0ab)
            k0_next = 0;

        if (!e1.has_on(pa) & (k0 != 1))
          if (e1ab)
            k0_next = 1;

        if (!e2.has_on(pa) & (k0 != 2))
          if (e2ab)
            k0_next = 2;

        if (k0_next == -1) {
          if (SHOW_QEDGE_DEBUG_INFO) {
            DEBUG_PRINT_INFO("  break, no k0_next");
          }
          break;
        }
        // do not trace over a non-constrained edge
        if (!_E_isNarrow(f0, k0_next)) {
          if (SHOW_QEDGE_DEBUG_INFO) {
            DEBUG_PRINT_INFO("  break, k0_next not constrained");
          }
          break;
        }

        // step into the adjacent face
        a = (complexi)_R_transfn(f0, k0_next) * a + (complexi)_T_integer_transfn(f0, k0_next);
        b = (complexi)_R_transfn(f0, k0_next) * b + (complexi)_T_integer_transfn(f0, k0_next);

        const int f1 = _TT(f0, k0_next);
        const int k1 = _TTi(f0, k0_next);

        if (SHOW_QEDGE_DEBUG_INFO) {
          DEBUG_PRINT_INFO("  f0= %6d : k0_next=%d : f1= %6d, k1=%d", //
                           f0, k0_next,                               //
                           // _F(f0, k0_next),                                        //
                           // _F(f0, (k0_next + 1) % 3),                              //
                           f1, k1);
        }

        if (f1 == -1) {
          if (SHOW_QEDGE_DEBUG_INFO) {
            DEBUG_PRINT_INFO("reached boundary.");
          }
          break;
        }

        if (!_F_isNarrow(f1)) {
          if (SHOW_QEDGE_DEBUG_INFO) {
            DEBUG_PRINT_INFO("reached white region.");
          }
          break;
        }
        // update the direction
        const int mm = _E_periodJumps(f0, k0_next);
        if (mm > 0) {
          r0 = (r0 + mm) % 4;
          if (SHOW_QEDGE_DEBUG_INFO) {
            DEBUG_PRINT_INFO("  update : r0 = %d", r0);
          }
        }
        f0 = f1;
        k0 = k1;

      } // end while

      if (SHOW_QEDGE_DEBUG_INFO)
        DEBUG_PRINT_INFO("point_is_inside_or_on_boundary=%d", point_is_inside_or_on_boundary);

      if (point_is_inside_or_on_boundary) {

        // find a qport with opposite direction
        const int r1_target = (r0 + 2) % 4;
        if (SHOW_QEDGE_DEBUG_INFO) {
          DEBUG_PRINT_INFO("  b      : u= %+d, v= %+d, r= %d", b.real(), b.imag(), r1_target);
        }

        // loop over all qports in f0
        for (auto qpi_next : _F_qports[f0]) {

          if (SHOW_QEDGE_DEBUG_INFO) {
            DEBUG_PRINT_INFO("  qp_next: u= %+d, v= %+d, r= %d    (%d)", //
                             _QP[qpi_next].u, _QP[qpi_next].v, _QP[qpi_next].pj, qpi_next);
          }

          const int r_next   = _QP[qpi_next].pj;
          const int qvi_next = _QP[qpi_next].qvi;

          if (r_next == r1_target            //
              && _QP[qpi_next].u == b.real() // direct comparison is ok since these are integers
              && _QP[qpi_next].v == b.imag()) {

            if (qvi_next == qvi) {
              DEBUG_PRINT_WARNING("Trying to connect q-ports (%d, %d) from the same q-vertex (%d)", qpi, qpi_next, qvi);
              continue;
            }

            if (QV_QV_connections(qvi, qvi_next)) {
              DEBUG_PRINT_WARNING("Q-vertices (%d, %d) already connected", qvi, qvi_next);
              if (edge_is_tangent)
                _QE[_QVPair_to_QE[std::make_pair(qvi, qvi_next)]].is_tangent = true;
              continue;
            }

            if (QP_connected(qpi_next)) {
              DEBUG_PRINT_WARNING("qp_next (%d) already connected, marking as irregular", qpi_next);
              _QP[qpi_next].is_irregular = true;
            }

            _QVPair_to_QE.insert(std::make_pair(std::make_pair(qvi, qvi_next), qei));
            _QVPair_to_QE.insert(std::make_pair(std::make_pair(qvi_next, qvi), qei));

            const Eigen::Vector2d edge_vector(_QV[qvi].x - _QV[qvi_next].x, //
                                              _QV[qvi].y - _QV[qvi_next].y);

            _QE.push_back(QEdge(qpi, qpi_next, qvi, qvi_next, edge_vector.norm(), //
                                edge_has_flipped_faces, edge_is_tangent));

            _QP[qpi].qei      = qei;
            _QP[qpi_next].qei = qei;

            _QP_QP_connections.push_back(std::make_pair(qpi, qpi_next));

            QP_connected(qpi)      = true;
            QP_connected(qpi_next) = true;

            QV_QV_connections(qvi, qvi_next) = true;
            QV_QV_connections(qvi_next, qvi) = true;

            qei++;
            break;
          }
        }
      } // if (point_is_inside_or_on_boundary)
    }   // for (all qports)

    DEBUG_PRINT_INFO("%zd q-edges", _QE.size());

    {
      // build QV-QE adjacency
      int qpi, qvi, qei;
      for (qvi = 0; qvi < n_qverts; ++qvi)
        _QV[qvi].qedges.clear();

      for (qpi = 0; qpi < n_qports; ++qpi) {
        qvi = _QP[qpi].qvi;
        qei = _QP[qpi].qei;
        if (qei > -1)
          _QV[qvi].qedges.push_back(qei);
      }
    }

    {
      // build QE-QE adjacency (only include tangent edges!)
      for (int qei = 0; qei < _QE.size(); ++qei)
        _QE[qei].qedges.clear();

      for (int qei0 = 0; qei0 < _QE.size(); ++qei0) {
        // loop over qedges adjacent to the first vertex
        for (auto qei1 : _QV[_QE[qei0].qvi0].qedges)
          if (qei0 != qei1 && _QE[qei1].is_tangent)
            _QE[qei0].qedges.push_back(qei1);

        // loop over qedges adjacent to the second vertex
        for (auto qei1 : _QV[_QE[qei0].qvi1].qedges)
          if (qei0 != qei1 && _QE[qei1].is_tangent)
            _QE[qei0].qedges.push_back(qei1);
      }
    }

    return true;
  }

} // namespace IGSV
