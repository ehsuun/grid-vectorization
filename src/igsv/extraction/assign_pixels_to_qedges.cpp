#include <igsv/extraction/assign_pixels_to_qedges.h>

#include <igsv/common/defs.h>
#include <igsv/extraction_entities/Pixel.h>
#include <igsv/extraction_entities/QEdge.h>
#include <igsv/extraction_entities/QPort.h>
#include <igsv/extraction_entities/QVertex.h>

#include <CGAL/Barycentric_coordinates_2/Segment_coordinates_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/nearest_neighbor_delaunay_2.h>
#include <CGAL/squared_distance_2.h>

#include <cmath>
#include <complex>

#define KNN 5 // number of neighbors to search

namespace IGSV {

  using complexd   = std::complex<double>;
  using EP_IC_K    = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Line_2     = EP_IC_K::Line_2;
  using Point_2    = EP_IC_K::Point_2;
  using SegCoords  = CGAL::Barycentric_coordinates::Segment_coordinates_2<EP_IC_K>;
  using uVb        = CGAL::Triangulation_vertex_base_with_info_2<unsigned, EP_IC_K>;
  using uFb        = CGAL::Triangulation_face_base_with_info_2<unsigned, EP_IC_K>;
  using IndexedTds = CGAL::Triangulation_data_structure_2<uVb, uFb>;
  using IndexedDT  = CGAL::Delaunay_triangulation_2<EP_IC_K, IndexedTds>;
  using VertexList = std::vector<IndexedDT::Vertex_handle>;

  // ==========================================================================

  void get_tangent_direction(const complexd& frame0,  //
                             const complexd& frame1,  //
                             const complexd& edgevec, //
                             complexd& tangent) {

    double ae = std::arg(edgevec);
    double a0 = std::arg(frame0);
    double a1 = std::arg(frame1);

    if (ae < 0)
      ae += 2 * M_PI;
    if (a0 < 0)
      a0 += 2 * M_PI;
    if (a1 < 0)
      a1 += 2 * M_PI;

    double d0 = std::min(std::abs(ae - a0), std::abs(ae - a0 - M_PI));
    while (d0 > 0.5 * M_PI)
      d0 = d0 - M_PI;
    d0 = std::abs(d0);

    double d1 = std::min(std::abs(ae - a1), std::abs(ae - a1 - M_PI));
    while (d1 > 0.5 * M_PI)
      d1 = d1 - M_PI;
    d1 = std::abs(d1);

    tangent = (d0 < d1) ? frame0 : frame1;
  }

  // ==========================================================================

  bool assign_pixels_to_qedges(const cv::Mat& _grey,                 //
                               const cv::Mat& _sw_01,                //
                               const Eigen::VectorXcd& _PX_XY,       //
                               const Eigen::VectorXi& _PX_I,         //
                               const Eigen::VectorXi& _PX_J,         //
                               const Eigen::VectorXi& _PX_fid,       //
                               const Eigen::MatrixXd& _X0_unit_comb, //
                               const Eigen::MatrixXd& _X1_unit_comb, //
                               const std::vector<QVertex>& _QV,      //
                               const std::vector<QPort>& _QP,        //
                               const std::vector<QEdge>& _QE,        //
                               const double _sw_avg,                 //
                               std::vector<Pixel>& _PX) {

    // prepare
    const int n_qedges      = _QE.size();
    const int n_dark_pixels = _PX_XY.rows();

    int n_samples = 0;

    _PX.resize(n_dark_pixels);

    // helper vars
    int qv0, qv1, f;
    double px, py, ax, ay, bx, by;
    complexd frame0, frame1, edgevec, tangent;
    VertexList knn_search;
    std::vector<double> bc, distances, tparams;
    std::vector<int> nearest_qedges;

    int pi, pj;
    unsigned char pixel_greyval_uint;
    double pixel_greyval_01;
    double pixel_w;

    // samples for finding nearest points (one per edge)
    std::vector<std::pair<Point_2, unsigned>> samples(n_qedges);

    for (int qei = 0; qei < n_qedges; ++qei) {
      qv0 = _QE[qei].qvi0;
      qv1 = _QE[qei].qvi1;

      samples[qei] = std::make_pair(Point_2(0.5 * (_QV[qv0].x + _QV[qv1].x), //
                                            0.5 * (_QV[qv0].y + _QV[qv1].y)),
                                    qei);
    }

    // triangulation over the samples
    IndexedDT tri;
    tri.insert(samples.begin(), samples.end());

    // loop over black pixels
    for (int k = 0; k < n_dark_pixels; k++) {

      // pixel position
      px = _PX_XY(k).real();
      py = _PX_XY(k).imag();
      pi = _PX_I(k);
      pj = _PX_J(k);

      // pixel weight
      pixel_w = (1.0 - (double)_grey.at<unsigned char>(pi, pj) / 255.) / _sw_01.at<double>(pi, pj);

      // pixel face
      f = _PX_fid(k);

      // pixel frame (unit)
      frame0 = complexd(_X0_unit_comb(f, 0), _X0_unit_comb(f, 1));
      frame1 = complexd(_X1_unit_comb(f, 0), _X1_unit_comb(f, 1));

      //// K nearest q-edges to this pixel
      knn_search.clear();
      CGAL::nearest_neighbors(tri, Point_2(px, py), KNN, std::back_inserter(knn_search));

      // collect edge indices
      nearest_qedges.clear();
      nearest_qedges.reserve(KNN);
      for (VertexList::const_iterator it = knn_search.begin(); it != knn_search.end(); it++)
        nearest_qedges.push_back((*it)->info());

      // compute the orthogonal distances
      distances.clear();
      distances.reserve(KNN);
      tparams.clear();
      tparams.reserve(KNN);

      for (int qei : nearest_qedges) {

        qv0 = _QE[qei].qvi0;
        qv1 = _QE[qei].qvi1;

        ax = _QV[qv0].x;
        ay = _QV[qv0].y;
        bx = _QV[qv1].x;
        by = _QV[qv1].y;

        // edge vector (unit)
        edgevec = complexd(bx - ax, by - ay) / _QE[qei].length;

        // determine which frame direction is tangent (= closer to edge vector)
        get_tangent_direction(frame0, frame1, edgevec, tangent);

        const Point_2 P(px, py); // pixel
        const Point_2 A(ax, ay); // startpoint
        const Point_2 B(bx, by); // endpoint
        const Line_2 line(A, B); // line AB
        SegCoords seg_coords(A, B);

        bc.clear();
        bc.reserve(2);
        seg_coords(line.projection(P), std::inserter(bc, bc.end()));
        double t = bc[1];

        // sanity check
        if (std::isnan(t)) {
          std::cerr << "    assign_pixels_to_qedges: t is nan! qei = " << qei << std::endl
                    << "                             P (" << px << ", " << py << ")" << std::endl
                    << "                             A (" << ax << ", " << ay << ")" << std::endl
                    << "                             B (" << bx << ", " << by << ")" << std::endl;
          return false;
        }

        // clip
        if (t < 0)
          t = 0;
        if (t > 1)
          t = 1;

        // compute the nearest point and the distance to it
        const Point_2 PP((1 - t) * ax + t * bx, (1 - t) * ay + t * by);

        distances.push_back(std::sqrt(CGAL::squared_distance(P, PP)));
        tparams.push_back(t);
      }

      // sort by distance
      std::vector<int> didx(distances.size());
      std::size_t n(0);
      std::generate(std::begin(didx), std::end(didx), [&] { return n++; });
      std::sort(std::begin(didx), std::end(didx), [&](int i1, int i2) { return distances[i1] < distances[i2]; });

      // store pixel samples
      int qei, qv0, qv1;
      double t, d;
      for (auto knn_i : didx) {

        if (distances[knn_i] > 4 * _sw_avg)
          continue; // ignore if too far

        qei = nearest_qedges[knn_i];
        t   = tparams[knn_i];
        d   = distances[knn_i];

        qv0 = _QE[qei].qvi0;
        qv1 = _QE[qei].qvi1;

        _PX[k].samples.push_back(Sample(k, qei, t, d,                           // pixel index
                                        px, py,                                 // pixel position
                                        tangent.real(), tangent.imag(),         // pixel tangent
                                        (1. - t) * _QV[qv0].x + t * _QV[qv1].x, // qedge footpoint, x
                                        (1. - t) * _QV[qv0].y + t * _QV[qv1].y, // qedge footpoint, y
                                        pixel_w                                 //
                                        ));
        n_samples++;

      } // loop over knn qedges
    }   // loop over pixels

    DEBUG_PRINT_INFO("%d samples", n_samples);

    return true;
  }

} // namespace IGSV
