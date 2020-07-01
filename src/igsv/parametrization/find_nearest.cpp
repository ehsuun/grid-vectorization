#include <igsv/parametrization/find_nearest.h>
#include <igsv/common/CornerType.h>

#include <CGAL/Barycentric_coordinates_2/Triangle_coordinates_2.h>

namespace IGSV {

  using EP_IC_K   = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Point_2   = EP_IC_K::Point_2;
  using TriCoords = CGAL::Barycentric_coordinates::Triangle_coordinates_2<EP_IC_K>;

  // ===========================================================================

  bool triangulation_cgal(const Eigen::MatrixXd& _V,                //
                          const std::vector<std::vector<int>>& _VT, //
                          IndexedDT& _tri) {

    //// construct
    std::vector<std::pair<EP_IC_K::Point_2, unsigned>> pp(_V.rows());
    for (int i = 0; i < _V.rows(); i++)
      pp[i] = std::make_pair(EP_IC_K::Point_2(_V(i, 0), _V(i, 1)), i);
    _tri.clear();
    _tri.insert(pp.begin(), pp.end());

    // sorted vertex-triangle indices
    std::vector<std::vector<int>> VT_sorted = _VT; // creates a copy
    for (auto&& vt : VT_sorted)
      std::sort(vt.begin(), vt.end()); // sort

    // for each face in tri, find the index of the corresponding row in F
    for (IndexedDT::Finite_faces_iterator it = _tri.finite_faces_begin(); it != _tri.finite_faces_end(); it++) {

      // get list of faces adjacent to each vertex
      const auto face_list_v0 = VT_sorted[it->vertex(0)->info()];
      const auto face_list_v1 = VT_sorted[it->vertex(1)->info()];
      const auto face_list_v2 = VT_sorted[it->vertex(2)->info()];

      // intersect lists for v0 and v1
      std::vector<int> face_list_v01;
      std::set_intersection(face_list_v0.begin(), face_list_v0.end(), face_list_v1.begin(), face_list_v1.end(),
                            back_inserter(face_list_v01));

      // intersect lists for v0, v1 and v2
      std::vector<int> face_list_v012;
      std::set_intersection(face_list_v2.begin(), face_list_v2.end(), face_list_v01.begin(), face_list_v01.end(),
                            back_inserter(face_list_v012));

      // there has to be a SINGLE face index left!
      if (face_list_v012.size() != 1) {
        std::cerr << "ERROR: Face index list -- there has to be only one! # = " //
                  << face_list_v012.size() << "("                               //
                  << it->vertex(0)->info() << ", "                              //
                  << it->vertex(1)->info() << ", "                              //
                  << it->vertex(2)->info() << ")" << std::endl;
        exit(EXIT_FAILURE);
      }
      assert(face_list_v012.size() == 1);

      // store it to the CGAL struct
      it->info() = face_list_v012[0];
    }

    return true;
  }

  // ===========================================================================

  int get_face_index(const IndexedDT& _tri, //
                     double px,             //
                     double py) {
    IndexedDT::Face_handle fh = _tri.locate(EP_IC_K::Point_2(px, py));
    return fh->info();
  }

  // ===========================================================================

  void get_barycentric_coords(const Eigen::MatrixXd& _V, //
                              const Eigen::MatrixXi& _F, //
                              int f,                     //
                              double px,                 //
                              double py,                 //
                              std::vector<double>& bc) {
    bc.clear();
    TriCoords tri_coords(Point_2(_V(_F(f, 0), 0), _V(_F(f, 0), 1)), //
                         Point_2(_V(_F(f, 1), 0), _V(_F(f, 1), 1)), //
                         Point_2(_V(_F(f, 2), 0), _V(_F(f, 2), 1)));
    tri_coords(Point_2(px, py), std::inserter(bc, bc.end()));
  }

  // ===========================================================================

  int get_nearest_vertex(const IndexedDT& _tri, //
                         double px,             //
                         double py,             //
                         double& vx,            //
                         double& vy) {
    IndexedDT::Vertex_handle vh = _tri.nearest_vertex(EP_IC_K::Point_2(px, py));
    vx                          = vh->point()[0];
    vy                          = vh->point()[1];
    return vh->info();
  }

  // ===========================================================================
  
  // std::complex<double> get_tangent(const IndexedDT& _tri,               //
  //                                  const std::complex<double>& p,       //
  //                                  const Eigen::MatrixXd& X0_unit_comb, //
  //                                  const Eigen::MatrixXd& X1_unit_comb, //
  //                                  const Eigen::MatrixXi& F_label) {
  //
  //   // which face?
  //   IndexedDT::Face_handle fh = _tri.locate(EP_IC_K::Point_2(p.real(), p.imag()));
  //
  //   const int f = fh->info();
  //   // face label?
  //   auto type = F_label(f);
  //   if (type == CORNER_TYPE_INTEGER_U)
  //     return std::complex<double>(X0_unit_comb(f, 0), X0_unit_comb(f, 1));
  //   else if (type == CORNER_TYPE_INTEGER_V)
  //     return std::complex<double>(X1_unit_comb(f, 0), X1_unit_comb(f, 1));
  //
  //   return std::complex<double>(0, 0);
  // }

  // ===========================================================================

  bool find_nearest(const Eigen::MatrixXd& _V,                //
                    const Eigen::MatrixXi& _F,                //
                    const std::vector<std::vector<int>>& _VT, //
                    const Eigen::VectorXcd& _PX_XY,           //
                    IndexedDT& _tri,                          //
                    Eigen::VectorXi& _PX_nearestID,           //
                    Eigen::VectorXcd& _PX_nearestXY,          //
                    Eigen::VectorXi& _PX_fid,                 //
                    Eigen::MatrixXd& _PX_bc) {

    triangulation_cgal(_V, _VT, _tri);

    int f;
    double px, py, vx, vy;
    std::vector<double> bc;

    const int n_dark = _PX_XY.rows();

    _PX_nearestID.resize(n_dark, 1);
    _PX_nearestXY.resize(n_dark, 1);
    _PX_fid.resize(n_dark, 1);
    _PX_bc.resize(n_dark, 3);

    for (unsigned k = 0; k < n_dark; k++) {
      px = _PX_XY(k).real();
      py = _PX_XY(k).imag();

      f = get_face_index(_tri, px, py);
      get_barycentric_coords(_V, _F, f, px, py, bc);

      _PX_fid(k) = f;
      _PX_bc.row(k) << bc[0], bc[1], bc[2];

      _PX_nearestID(k) = get_nearest_vertex(_tri, px, py, vx, vy);
      _PX_nearestXY(k) = std::complex<double>(vx, vy);
    }

    return true;
  }

  // ===========================================================================

} // namespace IGSV
