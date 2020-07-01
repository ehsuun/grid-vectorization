#include <igsv/parametrization/triangulate.h>
#include <igsv/common/defs.h>
#include <igsv/parametrization/pix_coords_from_xy.h>

#include <igl/barycenter.h>
#include <igl/triangle/triangulate.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/vertex_triangle_adjacency.h>

#include <iostream>
#include <set>
#include <sstream>

namespace IGSV {

  // ===========================================================================

  bool triangulate(const cv::Mat& _bw,                  //
                   const cv::Mat& _dt,                  //
                   double _maxarea,                     //
                   bool _adaptive,                      //
                   double _narrow_band_radius,          //
                   double _avg_stroke_width,            //
                   Eigen::MatrixXd& _V,                 //
                   Eigen::MatrixXi& _F,                 //
                   Eigen::MatrixXd& _UV_generated,      //
                   Eigen::MatrixXd& _BC,                //
                   Eigen::MatrixXi& _TT,                //
                   Eigen::MatrixXi& _TTi,               //
                   std::vector<std::vector<int>>& _VT,  //
                   std::vector<std::vector<int>>& _VTi, //
                   Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isBlack) {

    const int _w = _bw.cols;
    const int _h = _bw.rows;

    //// i) Dense triangulation
    // -- points (black pixels & four image corners)
    Eigen::MatrixXd P_corners;
    P_corners.resize(4, 2);
    P_corners << 0, 0, _w, 0, _w, _h, 0, _h;

    // -- edges (boundary)
    Eigen::MatrixXi E_corners;
    E_corners.resize(4, 2);
    E_corners << 0, 1, 1, 2, 2, 3, 3, 0;

    // -- hole (outside region)
    Eigen::MatrixXd H_outside;
    H_outside.resize(1, 2);
    H_outside << 10 * _w, 10 * _h;

    // -- generate
    std::stringstream dense_param;
    dense_param << "-Qq30a" << _maxarea;
    igl::triangle::triangulate(P_corners, E_corners, H_outside, dense_param.str(), _V, _F);

    //// iii) adaptive triangulation (if enabled)
    if (_adaptive) {

      // pixel coords
      double x, y;
      int i, j;

      std::set<int> vertex_list;

      const double base_scale = _narrow_band_radius * _avg_stroke_width;

      //// o) face-based
      // barycenters
      igl::barycenter(_V, _F, _BC);
      for (unsigned f = 0; f < _F.rows(); ++f) {
        x = _BC(f, 0);
        y = _BC(f, 1);
        pix_coords_from_xy(_w, _h, x, y, i, j);
        if ((double)_dt.at<unsigned char>(i, j) > base_scale)
          continue;
        vertex_list.insert(_F(f, 0));
        vertex_list.insert(_F(f, 1));
        vertex_list.insert(_F(f, 2));
      }

      const int n_dark = vertex_list.size();

      Eigen::MatrixXd P_tri;
      P_tri.resize(n_dark + 4, 2);
      int p = 0;
      for (auto v : vertex_list)
        P_tri.row(p++) << _V.row(v);
      P_tri.bottomRows(4) << P_corners;

      Eigen::MatrixXi E_tri = E_corners;
      E_tri.array() += n_dark; // shift, corners are at the end

      std::stringstream adaptive_param;
      adaptive_param << "-Qq30";
      igl::triangle::triangulate(P_tri, E_tri, H_outside, adaptive_param.str(), _V, _F);
    }

    //// IMPORTANT: V has only two cols!
    // append z=0 to avoid problems
    _V.conservativeResize(_V.rows(), 3);
    _V.col(2).fill(0);

    DEBUG_PRINT_INFO("nf = %zd", _F.rows());
    DEBUG_PRINT_INFO("nv = %zd", _V.rows());

    return process_mesh(_bw, _V, _F, _UV_generated, _BC, _TT, _TTi, _VT, _VTi, _F_isBlack);
  }

  // ===========================================================================

  bool process_mesh(const cv::Mat& _bw,                  //
                    const Eigen::MatrixXd& _V,           //
                    const Eigen::MatrixXi& _F,           //
                    Eigen::MatrixXd& _UV_generated,      //
                    Eigen::MatrixXd& _BC,                //
                    Eigen::MatrixXi& _TT,                //
                    Eigen::MatrixXi& _TTi,               //
                    std::vector<std::vector<int>>& _VT,  //
                    std::vector<std::vector<int>>& _VTi, //
                    Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isBlack) {

    const int _w = _bw.cols;
    const int _h = _bw.rows;

    //// generated uvs for texturing
    _UV_generated = _V.leftCols<2>();
    _UV_generated.col(0).array() /= _w;
    _UV_generated.col(1).array() /= _h;

    //// barycenters
    igl::barycenter(_V, _F, _BC);

    //// adjacency
    igl::triangle_triangle_adjacency(_F, _TT, _TTi);
    igl::vertex_triangle_adjacency(_V.rows(), _F, _VT, _VTi);

    //// which faces are black?
    _F_isBlack.resize(_F.rows(), 1);
    int i, j; // pix coords
    for (int f = 0; f < _F.rows(); f++) {
      pix_coords_from_xy(_w, _h, _BC(f, 0), _BC(f, 1), i, j);
      _F_isBlack(f) = _bw.at<unsigned char>(i, j) == 0;
    }

    DEBUG_PRINT_INFO("%zd black triangles", _F_isBlack.count());

    return true;
  }

  // ===========================================================================

} // namespace IGSV
