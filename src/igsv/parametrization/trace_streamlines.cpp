#include <igsv/common/CornerType.h>
#include <igsv/common/FaceData.h>
#include <igsv/parametrization/connect_triangles.h>
#include <igsv/parametrization/pix_coords_from_xy.h>
#include <igsv/parametrization/trace_streamlines.h>

#include <CoMISo/Utils/StopWatch.hh>
#include <boost/format.hpp>
#include <igsv/common/print_timings.h>
#include <iostream>

namespace IGSV {

  // ===========================================================================

  bool trace_streamlines(const Eigen::MatrixXi& _TT,                                //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isBlack,  //
                         const Eigen::MatrixXd& _F_barycenters,                     //
                         const Eigen::VectorXcd (&_F_frameField)[4],                //
                         const Eigen::MatrixXi& _E_periodJumps,                     //
                         const IndexedDT& _tri,                                     //
                         const cv::Mat& _bw,                                        //
                         double _sw_avg,                                            //
                         double _scale_multiplier,                                  //
                         double _tangent_ratio_threshold,                           //
                         double _tangent_ratio_streamlen,                           //
                         Eigen::VectorXi& _F_labels,                                //
                         std::vector<FaceData>& _FaceDataList) {

    // For a given frame, determine the tangent direction
    // by counting the number of black pixels along a line
    // in one of the directions of the frame field.
    // Tangent direction is the one which has greater number of black pixels.

    COMISO::StopWatch sw;

    double time__pointloc   = 0.;
    double time__connecttri = 0.;

    const int MIN_DILATION_SIZE = 1;

    long d_size = std::round(0.25 * _sw_avg);
    if (d_size < MIN_DILATION_SIZE)
      d_size = MIN_DILATION_SIZE;

    cv::Mat bw_streamlines;
    cv::erode(_bw, bw_streamlines,
              cv::getStructuringElement(cv::MORPH_ELLIPSE, cv::Size(2 * d_size + 1, 2 * d_size + 1),
                                        cv::Point(d_size, d_size)));

    const int numsamples = std::round(_tangent_ratio_streamlen * _sw_avg);

    const int nf       = _F_isBlack.rows();
    const int nf_black = _F_isBlack.count();

    const int h = bw_streamlines.rows;
    const int w = bw_streamlines.cols;

    _FaceDataList.clear();

    // barycenter position
    double bx, by;

    // pixel coords
    double pix_x, pix_y;
    int pix_i, pix_j;

    // current face handle
    IndexedDT::Face_handle fh;

    // current face index
    int f0;

    // current direction
    int mm0;

    // current position
    std::complex<double> p0;

    // current frame field vector
    std::complex<double> v0;
    bool update_v0;

    // next face index
    int f1;

    // next position
    std::complex<double> p1;

    // path between triangles
    std::vector<int> faceList;
    std::vector<int> edgeList;

    // tracing step
    const double dt = 1.0;

    // init a progress bar
    int n_progress_bars = 0;
    progress_bar(-1, nf, 10, n_progress_bars);

    int f_black = 0;

    for (int f = 0; f < nf; f++) {

      // skip white triangles
      if (!_F_isBlack(f))
        continue;

      // update progress bar
      progress_bar(f_black++, nf_black, 10, n_progress_bars);

      // barycenter position
      bx = _F_barycenters(f, 0);
      by = _F_barycenters(f, 1);

      // init face data
      FaceData fcdt(f, bx, by, numsamples);
      fh = IndexedDT::Face_handle();

      // for each of the four directions ...
      for (int mm_init = 0; mm_init < 4; ++mm_init) {

        // init
        p0        = std::complex<double>(bx, by);
        f0        = f;
        mm0       = mm_init;
        update_v0 = true;

        // stored positions of traced samples (numsamples x 3 matrix)
        fcdt.samples[mm_init].resize(fcdt.numsamples + 1, 3);
        fcdt.samples[mm_init].row(0) << bx, by, 0.0;

        for (int ns = 0; ns < fcdt.numsamples; ns++) {

          // 1. current frame field vector
          if (update_v0)
            v0 = _F_frameField[mm0](f0);

          // 2. step in the direction of the current frame field vector
          p1 = p0 + dt * v0;

          // 3. locate the new point in the triangulation
          sw.start();
          fh = _tri.locate(EP_IC_K::Point_2(p1.real(), p1.imag()), fh);
          time__pointloc += sw.stop() / 1000.0;

          // leaving the mesh, stop tracing
          if (_tri.is_infinite(fh)) {
            fcdt.samples[mm_init].conservativeResize(ns + 1, 3);
            break;
          }

          // index of the face in F
          f1 = fh->info();

          if (!_F_isNarrow(f1)) {
            fcdt.samples[mm_init].conservativeResize(ns + 1, 3);
            break;

            // leaving the narrow band, freeze the vector
            update_v0 = false;
          }

          // 4. compute the path from f0 to f1 >> faceList, edgeList = path(f0,f1)
          sw.start();
          // faceList and edgeList are cleared in connect_triangles
          connect_triangles(_TT, _tri, f0, f1, faceList, edgeList, p0.real(), p0.imag(), p1.real(), p1.imag(), fh);
          time__connecttri += sw.stop() / 1000.0;

          // 5. update
          // -- face
          f0 = f1;
          // -- position
          p0 = p1;
          // -- mmatch
          bool should_stop = false;
          for (int i = 0; i < edgeList.size(); i++) {
            const int _f = faceList[i];
            const int _k = edgeList[i];

            mm0 = (mm0 + _E_periodJumps(_f, _k)) % 4;
          }

          if (should_stop)
            break;

          // instead : store the path
          std::copy(faceList.begin(), faceList.end() - 1 /* exclude last face */,
                    std::back_inserter(fcdt.fpath[mm_init]));
          std::copy(edgeList.begin(), edgeList.end(), std::back_inserter(fcdt.epath[mm_init]));

          // 6. store
          pix_x = p1.real();
          pix_y = p1.imag();

          // -- position
          fcdt.samples[mm_init].row(ns + 1) << pix_x, pix_y, 0.0;

          // -- update black pixel count
          pix_coords_from_xy(w, h, pix_x, pix_y, pix_i, pix_j);

          if (bw_streamlines.at<unsigned char>(pix_i, pix_j) == 0)
            fcdt.c[mm_init % 2]++;
        }
      }

      // store face data
      _FaceDataList.push_back(fcdt);
    } // end loop over faces

    _F_labels.resize(nf, 1);
    _F_labels.fill(CORNER_TYPE_NO_INTEGERS);

    double ratio;

    for (auto const& fcdt : _FaceDataList) {

      ratio = (double)fcdt.c[0] / (double)(fcdt.c[0] + fcdt.c[1]);

      if (ratio < _tangent_ratio_threshold)
        _F_labels(fcdt.f) = CORNER_TYPE_INTEGER_V;

      else if (ratio > 1.0 - _tangent_ratio_threshold)
        _F_labels(fcdt.f) = CORNER_TYPE_INTEGER_U;
    }

    return true;
  }

} // namespace IGSV
