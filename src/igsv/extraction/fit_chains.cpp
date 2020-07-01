#include <igsv/extraction/fit_chains.h>

#include <igsv/common/defs.h>
#include <igsv/extraction_entities/Chain.h>
#include <igsv/extraction_entities/QEdge.h>
#include <igsv/extraction_entities/QVertex.h>

#include <alglib/src/interpolation.cpp>

#include <iostream>

namespace IGSV {

  using Matrix8d = Eigen::Matrix<double, 8, 8>;
  using Vector8d = Eigen::Matrix<double, 8, 1>;

  //===========================================================================

  int fit_cubic_spline_to_qverts(const std::vector<QVertex>& _QV,                         //
                                 const std::map<std::pair<int, int>, int>& _QVPair_to_QE, //
                                 const std::vector<QEdge>& _QE,                           //
                                 Chain& _chain,                                           //
                                 int cui,                                                 //
                                 int n_segments) {

    const int n_qverts = _chain.qverts.size();

    alglib::spline1dinterpolant spline_x, spline_y;

    //// Cubic spline, fit to qvertices
    {
      alglib::real_1d_array t, x, y, w; // samples = qvertices
      alglib::real_1d_array tc, xc, yc; // constraints
      alglib::integer_1d_array dc;      // types of constraints
      int n_constr;

      if (_chain.qvi0 == _chain.qvi1) {

        t.setlength(n_qverts + 1);
        x.setlength(n_qverts + 1);
        y.setlength(n_qverts + 1);
        w.setlength(n_qverts + 1);
        for (int i = 0; i < n_qverts; ++i) {
          t(i) = _chain.t[i];
          x(i) = _QV[_chain.qverts[i]].x;
          y(i) = _QV[_chain.qverts[i]].y;
          w(i) = 1.0; // uniform weights
        }

        t(n_qverts) = _chain.t[n_qverts];
        x(n_qverts) = _QV[_chain.qverts[0]].x;
        y(n_qverts) = _QV[_chain.qverts[0]].y;
        w(n_qverts) = 1.0;

        // n_constr = 2;
        // tc.setlength(n_constr);
        // xc.setlength(n_constr);
        // yc.setlength(n_constr);
        // dc.setlength(n_constr);
        //
        // tc(0) = 0.;
        // xc(0) = _QV[_chain.qverts[0]].x;
        // yc(0) = _QV[_chain.qverts[0]].y;
        // dc(0) = 0; // C0 constraint
        //
        // tc(1) = 1.;
        // xc(1) = _QV[_chain.qverts[0]].x;
        // yc(1) = _QV[_chain.qverts[0]].y;
        // dc(1) = 0; // C0 constraint

        alglib::spline1dbuildcatmullrom(t, x, n_qverts + 1, -1, 0.0, spline_x);
        alglib::spline1dbuildcatmullrom(t, y, n_qverts + 1, -1, 0.0, spline_y);
        // alglib::spline1dbuildcubic(t, x, n_qverts + 1, -1, 0., -1, 0., spline_x);
        // alglib::spline1dbuildcubic(t, y, n_qverts + 1, -1, 0., -1, 0., spline_y);
        _chain.fitting_cost = 0.;

      } else {

        t.setlength(n_qverts);
        x.setlength(n_qverts);
        y.setlength(n_qverts);
        w.setlength(n_qverts);
        for (int i = 0; i < n_qverts; ++i) {
          t(i) = _chain.t[i];
          x(i) = _QV[_chain.qverts[i]].x;
          y(i) = _QV[_chain.qverts[i]].y;
          w(i) = 1.0; // uniform weights
        }

        n_constr = 2;
        tc.setlength(n_constr);
        xc.setlength(n_constr);
        yc.setlength(n_constr);
        dc.setlength(n_constr);
        tc(0) = _chain.t[0];
        xc(0) = _QV[_chain.qverts[0]].x;
        yc(0) = _QV[_chain.qverts[0]].y;
        dc(0) = 0; // C0 constraint
        tc(1) = _chain.t[n_qverts - 1];
        xc(1) = _QV[_chain.qverts[n_qverts - 1]].x;
        yc(1) = _QV[_chain.qverts[n_qverts - 1]].y;
        dc(1) = 0; // C0 constraint

        // Rep -   report, same format as in LSFitLinearWC() subroutine.
        //             Following fields are set:
        //             * RMSError      rms error on the (X,Y).
        //             * AvgError      average error on the (X,Y).
        //             * AvgRelError   average relative error on the non-zero Y
        //             * MaxError      maximum error
        //                             NON-WEIGHTED ERRORS ARE CALCULATED
        alglib::spline1dfitreport rep_x, rep_y;
        alglib::ae_int_t info;
        alglib::spline1dfitcubicwc(t, x, w, n_qverts, tc, xc, dc, n_constr, n_segments + 3, info, spline_x, rep_x);
        alglib::spline1dfitcubicwc(t, y, w, n_qverts, tc, yc, dc, n_constr, n_segments + 3, info, spline_y, rep_y);
        // std::cout << "    x error: " << std::endl                         //
        //           << "    -- max    = " << rep_x.maxerror << std::endl    //
        //           << "    -- avg    = " << rep_x.avgerror << std::endl    //
        //           << "    -- rms    = " << rep_x.rmserror << std::endl    //
        //           << "    -- avgrel = " << rep_x.avgrelerror << std::endl //
        //     ;
        // std::cout << "    y error: " << std::endl                         //
        //           << "    -- max    = " << rep_y.maxerror << std::endl    //
        //           << "    -- avg    = " << rep_y.avgerror << std::endl    //
        //           << "    -- rms    = " << rep_y.rmserror << std::endl    //
        //           << "    -- avgrel = " << rep_y.avgrelerror << std::endl //
        ;
        _chain.fitting_cost = std::max(rep_x.maxerror, rep_y.maxerror);

        // alglib::spline1dbuildcatmullrom(t, x, spline_x);
        // alglib::spline1dbuildcatmullrom(t, y, spline_y);
        // alglib::spline1dbuildcubic(t, x, spline_x);
        // alglib::spline1dbuildcubic(t, y, spline_y);
        // _chain.fitting_cost = 0.;
      }
    }

    _chain.normalized_cost = _chain.fitting_cost;

    //// unpack
    {
      alglib::real_2d_array table_x, table_y;
      alglib::ae_int_t n;
      alglib::spline1dunpack(spline_x, n, table_x);
      alglib::spline1dunpack(spline_y, n, table_y);

      _chain.CubicSpline.resize(table_x.rows(), 10);
      for (int i = 0; i < table_x.rows(); ++i) {
        _chain.CubicSpline(i, 0) = table_x(i, 0);
        _chain.CubicSpline(i, 1) = table_x(i, 1);
        _chain.CubicSpline(i, 2) = table_x(i, 2);
        _chain.CubicSpline(i, 3) = table_x(i, 3);
        _chain.CubicSpline(i, 4) = table_x(i, 4);
        _chain.CubicSpline(i, 5) = table_x(i, 5);
        _chain.CubicSpline(i, 6) = table_y(i, 2);
        _chain.CubicSpline(i, 7) = table_y(i, 3);
        _chain.CubicSpline(i, 8) = table_y(i, 4);
        _chain.CubicSpline(i, 9) = table_y(i, 5);
      }
    }

    return _chain.CubicSpline.rows(); // return the number of segments
  }

  //===========================================================================

  void coeffs_monomial_to_bernstein(Eigen::Vector4d& _c, double a, double b) {

    // in-place conversion to Bernstein basis from the monomial basis

    const double b_a = b - a;

    // 1. convert to monomial basis with normalized parameter t in [0,1]
    _c(0) = _c(0) + _c(1) * a + _c(2) * a * a + _c(3) * a * a * a;
    _c(1) = b_a * (_c(1) + 2 * _c(2) * a + 3 * _c(3) * a * a);
    _c(2) = b_a * b_a * (_c(2) + 3 * _c(3) * a);
    _c(3) = b_a * b_a * b_a * _c(3);

    // 2. convert to Bernstein basis, t in [0,1]
    _c(3) = _c(0) + _c(1) + _c(2) + _c(3);
    _c(2) = _c(0) + _c(1) * 2. / 3. + _c(2) / 3.;
    _c(1) = _c(0) + _c(1) / 3.;
  }

  //===========================================================================

  bool convert_splines_to_bezier_form(std::vector<Chain>& _CH) {

    for (auto& _chain : _CH) {

      _chain.CP_spline.resize(4, _chain.CubicSpline.rows());

      for (int i = 0; i < _chain.CubicSpline.rows(); ++i) {

        Eigen::Vector4d _cx, _cy;
        for (int k = 0; k < 4; ++k) {
          _cx(k) = _chain.CubicSpline(i, 2 + k);
          _cy(k) = _chain.CubicSpline(i, 6 + k);
        }

        const double a = _chain.CubicSpline(i, 0);
        const double b = _chain.CubicSpline(i, 1);

        coeffs_monomial_to_bernstein(_cx, 0.0, b - a);
        coeffs_monomial_to_bernstein(_cy, 0.0, b - a);

        for (int k = 0; k < 4; ++k)
          _chain.CP_spline(k, i) = std::complex<double>(_cx(k), _cy(k));
      }
    }

    return true;
  }

  //===========================================================================

  bool fit_chains(const std::vector<QVertex>& _QV,                         //
                  const std::map<std::pair<int, int>, int>& _QVPair_to_QE, //
                  std::vector<QEdge>& _QE,                                 //
                  std::vector<Chain>& _CH,                                 //
                  double _avg_stroke_width,                                //
                  double _narrow_band_radius) {

    const double max_fitting_cost = _avg_stroke_width * _narrow_band_radius;

    const int max_fitting_iters = 10; // 1 means do not split

#define HISTOGRAM_SIZE 4

    Eigen::VectorXi histogram;
    histogram.resize(HISTOGRAM_SIZE, 1);

    int cui = 0;
    for (auto& _chain : _CH) {

      const int n_qverts = _chain.qverts.size();

      // check coverage
      histogram.fill(0);
      for (const auto& sample : _chain.samples)
        for (int k = 1; k <= HISTOGRAM_SIZE; k++)
          if (sample.t_chain < (double)k / (double)HISTOGRAM_SIZE || k == HISTOGRAM_SIZE) {
            histogram[k - 1]++;
            break;
          }

      if ((histogram.array() == 0).any())
        continue; // skip this chain

      try {
        int n_segments = 1;
        int n_iters    = 0;

        while (_chain.normalized_cost > max_fitting_cost && n_segments < n_qverts && n_iters++ < max_fitting_iters) {
          n_segments = fit_cubic_spline_to_qverts(_QV, _QVPair_to_QE, _QE, _chain, cui, n_segments);
          n_segments++; // try to increase the number of segments
        }

        convert_splines_to_bezier_form(_CH);

      } catch (alglib::ap_error& e) {
        DEBUG_PRINT_WARNING("curve #%d: alglib exception caught", cui);
      }

      cui++;
    }

    return true;
  }

} // namespace IGSV
