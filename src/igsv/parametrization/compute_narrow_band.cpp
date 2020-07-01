#include <igsv/parametrization/pix_coords_from_xy.h>
#include <Eigen/Core>
#include <opencv2/opencv.hpp>

namespace IGSV {

  bool compute_narrow_band(const Eigen::MatrixXd& _V,                                   //
                           const Eigen::MatrixXi& _F,                                   //
                           const Eigen::MatrixXi& _TT,                                  //
                           const Eigen::MatrixXi& _TTi,                                 //
                           const std::vector<std::vector<int>>& _VT,                    //
                           const std::vector<std::vector<int>>& _VTi,                   //
                           const Eigen::MatrixXd& _BC,                                  //
                           const cv::Mat& _bw,                                          //
                           const cv::Mat& _dt,                                          //
                           const cv::Mat& _detail_mask,                                 //
                           double _mask_factor,                                         //
                           double _narrow_band_radius,                                  //
                           double _avg_stroke_width,                                    //
                           Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isBlack,          //
                           Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrow,         //
                           Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrowBoundary, //
                           Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow,         //
                           Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrow,         //
                           Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrowBoundary) {

    const int w = _bw.cols;
    const int h = _bw.rows;

    int i, j;

    // determine black vertices
    _V_isBlack.resize(_V.rows(), 1);
    _V_isBlack.fill(false);
    for (int v = 0; v < _V.rows(); ++v) {
      pix_coords_from_xy(w, h, _V(v, 0), _V(v, 1), i, j);
      _V_isBlack(v) = _bw.at<unsigned char>(i, j) == 0;
    }

    // determine narrow band faces and vertices
    _F_isNarrow.resize(_F.rows(), 1);
    _F_isNarrow.fill(false);
    _V_isNarrow.resize(_V.rows(), 1);
    _V_isNarrow.fill(false);

    const double distance_threshold = _narrow_band_radius * _avg_stroke_width;

    double mask_multiplier;

    for (int f = 0; f < _F.rows(); ++f) {
      pix_coords_from_xy(w, h, _BC(f, 0), _BC(f, 1), i, j);

      mask_multiplier = 1.;
      if (_detail_mask.at<unsigned char>(i, j) == 255) {
        // coarser scale, wider band : multiply the threshold (bigger threshold)
        mask_multiplier = _mask_factor;
      } else if (_detail_mask.at<unsigned char>(i, j) == 0) {
        // finer scale, narrower band : divide the threshold (smaller threshold)
        mask_multiplier = 1. / _mask_factor;
      }

      if ((double)_dt.at<unsigned char>(i, j) > distance_threshold * mask_multiplier)
        continue;
      _F_isNarrow(f) = true;
    }

    // remove isolated triangles (to prevent parametrization degeneracies)
    for (int f = 0; f < _F.rows(); ++f)
      if (_F_isNarrow(f)) {
        // is at least one neighboring triangle in the narrow band?
        bool has_adjacent_triangles_in_nb = false;
        for (int k = 0; k < 3; ++k)
          if (_F_isNarrow(_TT(f, k))) {
            has_adjacent_triangles_in_nb = true;
            break;
          }
        _F_isNarrow(f) = has_adjacent_triangles_in_nb;
      }

    // vertices in the narrow band
    for (int f = 0; f < _F.rows(); ++f)
      if (_F_isNarrow(f)) {
        _V_isNarrow(_F(f, 0)) = true;
        _V_isNarrow(_F(f, 1)) = true;
        _V_isNarrow(_F(f, 2)) = true;
      }

    //// repair the mesh : fill small regions and remove non-manifold vertices
    int n_filled;
    int n_nonmanifold;
    int n_repair_iter = 0;
    do {

      // fill holes : face is in the narrow band if all its vertices are in the narrow band
      n_filled = 0;
      for (int f = 0; f < _F.rows(); ++f)
        if (!_F_isNarrow(f))
          if (_V_isNarrow(_F(f, 0)) & _V_isNarrow(_F(f, 1)) & _V_isNarrow(_F(f, 2))) {
            _F_isNarrow(f) = true;
            n_filled++;
          }

      // repair non-manifold vertices
      n_nonmanifold = 0;
      for (int v = 0; v < _V.rows(); ++v) {

        const int f_start = _VT[v][0];
        const int k_start = _VTi[v][0];

        int f0 = f_start;
        int k0 = k_start;
        int f1, k1;

        int n_changes = 0;
        bool current;
        bool first = true;

        do {
          f1 = _TT(f0, k0);
          k1 = (_TTi(f0, k0) + 1) % 3;
          f0 = f1;
          k0 = k1;

          if (f0 == -1)
            break;

          if (first) {
            current = _F_isNarrow(f0);
            first   = false;
          } else {
            if (current != _F_isNarrow(f0)) {
              n_changes++;
              current = _F_isNarrow(f0);
            }
          }

        } while (f0 != f_start);

        if (n_changes > 2) {
          n_nonmanifold++;
          for (auto f : _VT[v])
            _F_isNarrow(f) = true;
        }
      }

      // std::cout << "    repair the mesh: n_filled = " << n_filled << ", n_nonmanifold = " << n_nonmanifold << std::endl;

    } while ((n_filled > 0 || n_nonmanifold > 0) && ++n_repair_iter < 10); // end repair mesh

    // determine narrow band edges
    _E_isNarrow.resize(_F.rows(), 3);
    _E_isNarrow.fill(false);

    _V_isNarrowBoundary.resize(_V.rows(), 1);
    _V_isNarrowBoundary.fill(false);

    _E_isNarrowBoundary.resize(_F.rows(), 3);
    _E_isNarrowBoundary.fill(false);

    Eigen::Matrix<bool, Eigen::Dynamic, 3> marked;
    marked.resize(_F.rows(), 3);
    marked.fill(false);

    int f0, k0, f1, k1;

    for (f0 = 0; f0 < _F.rows(); ++f0)
      for (k0 = 0; k0 < 3; ++k0) {

        if (marked(f0, k0))
          continue; // skip if already processed

        marked(f0, k0) = true; // mark as done

        f1 = _TT(f0, k0);
        k1 = _TTi(f0, k0);

        if (f1 == -1)
          continue; // skip if boundary

        marked(f1, k1) = true; // mark as done

        if (_F_isNarrow(f0) & _F_isNarrow(f1)) {
          _E_isNarrow(f0, k0) = true;
          _E_isNarrow(f1, k1) = true;
        }

        // boundary : one face inside, the other outside
        if ((_F_isNarrow(f0) & !_F_isNarrow(f1)) | (!_F_isNarrow(f0) & _F_isNarrow(f1))) {
          _V_isNarrowBoundary(_F(f0, k0))           = true;
          _V_isNarrowBoundary(_F(f0, (k0 + 1) % 3)) = true;
          _E_isNarrowBoundary(f0, k0)               = true;
          _E_isNarrowBoundary(f1, k1)               = true;
        }
      }

    return true;
  }

} // namespace IGSV
