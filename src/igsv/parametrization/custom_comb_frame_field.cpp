#include <Eigen/Core>
#include <deque>
#include <igl/PI.h>
#include <iostream>
#include <utility> // std::pair
#include <vector>

double Sign(double a) { return (double)((a > 0) ? +1 : -1); }

namespace IGSV {

  // ===========================================================================

  void convert_to_vectors(const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                          const Eigen::VectorXd& _X0_arg,                            //
                          const Eigen::VectorXd& _X1_arg,                            //
                          const Eigen::VectorXd& _X0_mag,                            //
                          const Eigen::VectorXd& _X1_mag,                            //
                          Eigen::MatrixXd& _X0_unit,                                 //
                          Eigen::MatrixXd& _X1_unit,                                 //
                          Eigen::MatrixXd& _X0_nonunit,                              //
                          Eigen::MatrixXd& _X1_nonunit) {

    //// Store both unit and non-unit frame field vectors
    _X0_unit.setZero(_X0_arg.rows(), 3);
    _X1_unit.setZero(_X0_arg.rows(), 3);
    _X0_nonunit.setZero(_X0_arg.rows(), 3);
    _X1_nonunit.setZero(_X0_arg.rows(), 3);

    double ux, uy, vx, vy;
    for (int f = 0; f < _X0_arg.rows(); ++f)
      if (_F_isNarrow(f)) {
        ux = std::cos(_X0_arg(f));
        uy = std::sin(_X0_arg(f));
        vx = std::cos(_X1_arg(f));
        vy = std::sin(_X1_arg(f));
        _X0_unit.row(f) << ux, uy, 0.0;
        _X1_unit.row(f) << vx, vy, 0.0;
        _X0_nonunit.row(f) << _X0_mag(f) * ux, _X0_mag(f) * uy, 0.0;
        _X1_nonunit.row(f) << _X1_mag(f) * vx, _X1_mag(f) * vy, 0.0;
      }
  }

  // ===========================================================================

  void bisector_field(const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                      const Eigen::MatrixXd& _X0_unit,                           //
                      const Eigen::MatrixXd& _X1_unit,                           //
                      Eigen::MatrixXd& _BIS0_unit,                               //
                      Eigen::MatrixXd& _BIS1_unit) {

    // adapted from igl::compute_frame_field_bisectors

    _BIS0_unit.setZero(_F_isNarrow.rows(), 3);
    _BIS1_unit.setZero(_F_isNarrow.rows(), 3);
    for (unsigned i = 0; i < _F_isNarrow.rows(); ++i)
      if (_F_isNarrow(i)) {
        // project onto the tangent plane and convert to angle
        // Convert to angle
        double a0 = atan2(_X0_unit(i, 1), _X0_unit(i, 0));
        // make it positive by adding some multiple of 2pi
        a0 += std::ceil(std::max(0., -a0) / (igl::PI * 2.)) * (igl::PI * 2.);
        // take modulo 2pi
        a0        = fmod(a0, (igl::PI * 2.));
        double a1 = atan2(_X1_unit(i, 1), _X1_unit(i, 0));
        // make it positive by adding some multiple of 2pi
        a1 += std::ceil(std::max(0., -a1) / (igl::PI * 2.)) * (igl::PI * 2.);
        // take modulo 2pi
        a1 = fmod(a1, (igl::PI * 2.));

        double b0 = (a0 + a1) / 2.0;
        // make it positive by adding some multiple of 2pi
        b0 += std::ceil(std::max(0., -b0) / (igl::PI * 2.)) * (igl::PI * 2.);
        // take modulo 2pi
        b0 = fmod(b0, (igl::PI * 2.));

        double b1 = b0 + (igl::PI / 2.);
        // make it positive by adding some multiple of 2pi
        b1 += std::ceil(std::max(0., -b1) / (igl::PI * 2.)) * (igl::PI * 2.);
        // take modulo 2pi
        b1 = fmod(b1, (igl::PI * 2.));

        _BIS0_unit.row(i) << cos(b0), sin(b0), 0.0;
        _BIS1_unit.row(i) << cos(b1), sin(b1), 0.0;
      }
  }

  // ===========================================================================

  void comb_bisectors(const Eigen::MatrixXi& _TT,                                //
                      const Eigen::MatrixXi& _TTi,                               //
                      const Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrow, //
                      const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                      const Eigen::MatrixXd& _BIS0_unit,                         //
                      const Eigen::MatrixXd& _BIS1_unit,                         //
                      Eigen::MatrixXd& _BIS0_unit_comb,                          //
                      Eigen::MatrixXd& _BIS1_unit_comb,                          //
                      Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_customCut,      //
                      std::vector<std::pair<int, int>>& _TT_dualConnections) {

    // Comb the bisector field
    _BIS0_unit_comb.setZero(_F_isNarrow.rows(), 3);
    _BIS1_unit_comb.setZero(_F_isNarrow.rows(), 3);
    _BIS0_unit_comb << _BIS0_unit;
    _BIS1_unit_comb << _BIS1_unit;

    // mark all edges as cut edges
    _E_customCut.resize(_F_isNarrow.rows(), 3);
    _E_customCut.fill(true);
    _TT_dualConnections.clear();

    // only process the narrow band
    Eigen::Matrix<bool, Eigen::Dynamic, 1> f_visited;
    f_visited.resize(_F_isNarrow.rows(), 1);
    for (int f = 0; f < _F_isNarrow.rows(); ++f)
      f_visited(f) = !_F_isNarrow(f);

    std::deque<int> queue;

    while (!f_visited.all()) {

      // there are faces that were not yet visited

      // find such face
      int f_start = 0;
      while (f_visited(f_start))
        f_start++;

      // add it as a seed
      queue.push_back(f_start);
      f_visited(f_start) = true;

      while (!queue.empty()) {

        const int f0 = queue.front();
        queue.pop_front();

        for (int k0 = 0; k0 < 3; k0++) {

          if (!_E_isNarrow(f0, k0))
            continue; // leaving the narrow band, stop

          const int f1 = _TT(f0, k0);  // the other face
          const int k1 = _TTi(f0, k0); // the other face

          if (f1 == -1)
            continue; // skip if boundary
          if (f_visited(f1))
            continue; // skip if already visited

          _E_customCut(f0, k0) = false;
          _E_customCut(f1, k1) = false;
          _TT_dualConnections.push_back(std::make_pair(f0, f1));

          // which of the two is closer?
          const double d0 = _BIS0_unit_comb.row(f1).dot(_BIS0_unit_comb.row(f0));
          const double d1 = _BIS1_unit_comb.row(f1).dot(_BIS0_unit_comb.row(f0));

          // store the closer one as the first vector
          if (fabs(d0) >= fabs(d1))
            _BIS0_unit_comb.row(f1) << _BIS0_unit_comb.row(f1) * Sign(d0);
          else
            _BIS0_unit_comb.row(f1) << _BIS1_unit_comb.row(f1) * Sign(d1);

          // the second vector is obtained by rotating the first vector
          _BIS1_unit_comb.row(f1) << -_BIS0_unit_comb(f1, 1), _BIS0_unit_comb(f1, 0), 0.0;

          f_visited(f1) = true;
          queue.push_back(f1);
        }
      }
    }

    for (int f = 0; f < _F_isNarrow.rows(); ++f)
      if (!f_visited(f))
        std::cerr << "WARNING: face " << f << " was *not* visited!\n";

    // flip to get positive CGAL orientation
    _BIS1_unit_comb.array() *= -1;
  }

  // ===========================================================================

  void period_jumps(const Eigen::MatrixXi& _TT,                                //
                    const Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrow, //
                    const Eigen::MatrixXd& _BIS0_unit_comb,                    //
                    const Eigen::MatrixXd& _BIS1_unit_comb,                    //
                    Eigen::MatrixXi& _E_periodJumps) {

    // adapted from igl::cross_field_mismatch

    _E_periodJumps.setZero(_E_isNarrow.rows(), 3);
    int f1, pjump;
    double angle_diff;
    const double pi_2 = igl::PI / 2.0;

    for (int f0 = 0; f0 < _E_isNarrow.rows(); ++f0)
      for (int k0 = 0; k0 < 3; ++k0)
        if (_E_isNarrow(f0, k0)) {

          f1 = _TT(f0, k0);
          if (f1 == -1)
            continue;

          angle_diff = atan2(_BIS0_unit_comb.row(f1).dot(_BIS1_unit_comb.row(f0)),
                             _BIS0_unit_comb.row(f1).dot(_BIS0_unit_comb.row(f0)));

          pjump = (int)std::floor((angle_diff / pi_2) + 0.5);
          if (pjump >= 0)
            pjump = pjump % 4;
          else
            pjump = (-(3 * pjump)) % 4;

          _E_periodJumps(f0, k0) = pjump;
        }
  }

  // ===========================================================================

  void store_singularities(const Eigen::MatrixXi& _F,                                 //
                           const std::vector<std::vector<int>>& _VT,                  //
                           const std::vector<std::vector<int>>& _VTi,                 //
                           const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrow, //
                           const Eigen::MatrixXi& _E_periodJumps,                     //
                           Eigen::VectorXi& _V_singIndex) {

    _V_singIndex.setZero(_V_isNarrow.rows(), 1);

    for (int v = 0; v < _V_isNarrow.rows(); v++)
      if (_V_isNarrow[v]) {
        int mismatch = 0;
        for (unsigned int i = 0; i < _VT[v].size(); i++) {
          // look for the vertex
          int j = -1;
          for (unsigned z = 0; z < 3; ++z)
            if (_F(_VT[v][i], z) == v)
              j = z;
          assert(j != -1);
          mismatch += _E_periodJumps(_VT[v][i], j);
        }
        mismatch        = mismatch % 4;
        _V_singIndex(v) = mismatch;
      }
  }

  // ===========================================================================

  void comb_from_bisectors(const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                           const Eigen::MatrixXd& _X0,                                //
                           const Eigen::MatrixXd& _X1,                                //
                           const Eigen::MatrixXd& _BIS0_unit_comb,                    //
                           const Eigen::MatrixXd& _BIS1_unit_comb,                    //
                           Eigen::MatrixXd& _X0_comb,                                 //
                           Eigen::MatrixXd& _X1_comb) {

    //// adapted from igl::comb_frame_field

    _X0_comb.setZero(_BIS0_unit_comb.rows(), 3);
    _X1_comb.setZero(_BIS1_unit_comb.rows(), 3);

    for (unsigned i = 0; i < _X0.rows(); ++i)
      if (_F_isNarrow(i)) {

        Eigen::Matrix<double, 4, 3> DIRs;
        DIRs << _X0.row(i), -_X0.row(i), _X1.row(i), -_X1.row(i);

        std::vector<double> a(4);

        double a_combed = atan2(_BIS0_unit_comb(i, 1), _BIS0_unit_comb(i, 0));

        // center on the combed sector center
        for (unsigned j = 0; j < 4; ++j) {
          a[j] = atan2(DIRs(j, 1), DIRs(j, 0)) - a_combed;
          // make it positive by adding some multiple of 2pi
          a[j] += std::ceil(std::max(0., -a[j]) / (igl::PI * 2.)) * (igl::PI * 2.);
          // take modulo 2pi
          a[j] = fmod(a[j], (igl::PI * 2.));
        }
        // now the max is u and the min is v

        int m = std::min_element(a.begin(), a.end()) - a.begin();
        int M = std::max_element(a.begin(), a.end()) - a.begin();

        // _F_isCombingOk(i) = ((m >= 0 && m <= 1) && (M >= 2 && M <= 3)) || ((m >= 2 && m <= 3) && (M >= 0 && M <= 1));

        _X0_comb.row(i) = DIRs.row(m);
        _X1_comb.row(i) = DIRs.row(M);
      }

    // flip to get positive CGAL orientation
    _X1_comb.array() *= -1;
  }

  // ===========================================================================

  bool custom_comb_frame_field(const Eigen::MatrixXd& _V,                                 //
                               const Eigen::MatrixXi& _F,                                 //
                               const Eigen::MatrixXi& _TT,                                //
                               const Eigen::MatrixXi& _TTi,                               //
                               const std::vector<std::vector<int>>& _VT,                  //
                               const std::vector<std::vector<int>>& _VTi,                 //
                               const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _V_isNarrow, //
                               const Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_isNarrow, //
                               const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                               const Eigen::VectorXd& _X0_arg,                            //
                               const Eigen::VectorXd& _X1_arg,                            //
                               const Eigen::VectorXd& _X0_mag,                            //
                               const Eigen::VectorXd& _X1_mag,                            //
                               Eigen::VectorXi& _V_singIndex,                             //
                               Eigen::MatrixXi& _E_periodJumps,                           //
                               Eigen::Matrix<bool, Eigen::Dynamic, 3>& _E_customCut,      //
                               std::vector<std::pair<int, int>>& _TT_dualConnections,     //
                               Eigen::MatrixXd& _BIS0_unit,                               //
                               Eigen::MatrixXd& _BIS1_unit,                               //
                               Eigen::MatrixXd& _BIS0_unit_comb,                          //
                               Eigen::MatrixXd& _BIS1_unit_comb,                          //
                               Eigen::MatrixXd& _X0_unit,                                 //
                               Eigen::MatrixXd& _X1_unit,                                 //
                               Eigen::MatrixXd& _X0_nonunit,                              //
                               Eigen::MatrixXd& _X1_nonunit,                              //
                               Eigen::MatrixXd& _X0_unit_comb,                            //
                               Eigen::MatrixXd& _X1_unit_comb,                            //
                               Eigen::MatrixXd& _X0_nonunit_comb,                         //
                               Eigen::MatrixXd& _X1_nonunit_comb,                         //
                               Eigen::VectorXcd (&_F_frameField)[4]) {

    convert_to_vectors(_F_isNarrow, // input
                       _X0_arg,     //
                       _X1_arg,     //
                       _X0_mag,     //
                       _X1_mag,     //
                       _X0_unit,    // output
                       _X1_unit,    //
                       _X0_nonunit, //
                       _X1_nonunit);

    bisector_field(_F_isNarrow, // input
                   _X0_unit,    //
                   _X1_unit,    //
                   _BIS0_unit,  // output
                   _BIS1_unit);

    comb_bisectors(_TT,             // input
                   _TTi,            //
                   _E_isNarrow,     //
                   _F_isNarrow,     //
                   _BIS0_unit,      //
                   _BIS1_unit,      //
                   _BIS0_unit_comb, // output
                   _BIS1_unit_comb, //
                   _E_customCut,    //
                   _TT_dualConnections);

    period_jumps(_TT,             // input
                 _E_isNarrow,     //
                 _BIS0_unit_comb, //
                 _BIS1_unit_comb, //
                 _E_periodJumps   // output
    );

    store_singularities(_F,             // input
                        _VT,            //
                        _VTi,           //
                        _V_isNarrow,    //
                        _E_periodJumps, //
                        _V_singIndex    // output
    );

    comb_from_bisectors(_F_isNarrow,     // input
                        _X0_unit,        //
                        _X1_unit,        //
                        _BIS0_unit_comb, //
                        _BIS1_unit_comb, //
                        _X0_unit_comb,   // output
                        _X1_unit_comb);

    comb_from_bisectors(_F_isNarrow,      // input
                        _X0_nonunit,      //
                        _X1_nonunit,      //
                        _BIS0_unit_comb,  //
                        _BIS1_unit_comb,  //
                        _X0_nonunit_comb, // output
                        _X1_nonunit_comb);

    for (int k = 0; k < 4; k++)
      _F_frameField[k].resize(_F.rows(), 1);

    for (int f = 0; f < _F.rows(); f++) {
      _F_frameField[0](f) = std::complex<double>(+_X0_unit_comb(f, 0), +_X0_unit_comb(f, 1));
      _F_frameField[1](f) = std::complex<double>(+_X1_unit_comb(f, 0), +_X1_unit_comb(f, 1));
      _F_frameField[2](f) = std::complex<double>(-_X0_unit_comb(f, 0), -_X0_unit_comb(f, 1));
      _F_frameField[3](f) = std::complex<double>(-_X1_unit_comb(f, 0), -_X1_unit_comb(f, 1));
    }

    return true;
  }

} // namespace IGSV
