#include <Eigen/Core>

#include <iostream>

#include <igsv/parametrization/process_uv.h>

namespace IGSV {

  // ===========================================================================

  bool sanitize_parametrization(Eigen::MatrixXd& _UV) {

    double u, v;
    for (int i = 0; i < _UV.rows(); i++) {

      u = _UV(i, 0);
      v = _UV(i, 1);

      if (std::abs(u - std::round(u)) < 1e-8)
        _UV(i, 0) = std::round(u);

      if (std::abs(v - std::round(v)) < 1e-8)
        _UV(i, 1) = std::round(v);
    }

    return true;
  }

  // ===========================================================================

  bool precompute_transition_fn(const Eigen::MatrixXd& _UV,            //
                                const Eigen::MatrixXi& _FUV,           //
                                const Eigen::MatrixXi& _TT,            //
                                const Eigen::MatrixXi& _TTi,           //
                                const Eigen::MatrixXi& _E_periodJumps, //
                                Eigen::MatrixXcd& _R_transfn,          //
                                Eigen::MatrixXcd& _T_transfn,          //
                                Eigen::MatrixXcd& _T_integer_transfn,  //
                                Eigen::MatrixXcd& _T_decimal_transfn,  //
                                Eigen::Matrix<bool, Eigen::Dynamic, 3>& _flip_transfn) {

    const int nf = _TT.rows();

    _R_transfn.resize(nf, 3);
    _T_transfn.resize(nf, 3);
    _T_integer_transfn.resize(nf, 3);
    _T_decimal_transfn.resize(nf, 3);
    _flip_transfn.resize(nf, 3);

    for (int f0 = 0; f0 < nf; f0++) {
      for (int k0 = 0; k0 < 3; k0++) {

        // adjacent face
        const int f1 = _TT(f0, k0);
        const int k1 = _TTi(f0, k0);

        if (f1 == -1)
          continue;

        // rotational part
        const int pjump = _E_periodJumps(f0, k0);
        std::complex<double> r01;
        switch (pjump) {
        case 0:
          r01 = std::complex<double>(+1, 0);
          break;
        case 1:
          r01 = std::complex<double>(0, +1);
          break;
        case 2:
          r01 = std::complex<double>(-1, 0);
          break;
        default:
          r01 = std::complex<double>(0, -1);
        }

        //// indices of vertices from both sides of the edge
        // std::cout << "f0=" << f0 << ",k0=" << k0 << ",f1=" << f1 << ",k1=" << k1 //
        //           << " | " << _FUV.rows() << "," << _FUV.cols() << std::endl;
        // std::cout << "    ip0" << std::endl;
        const int ip0 = _FUV(f0, k0);
        // std::cout << "    iq0" << std::endl;
        const int iq0 = _FUV(f0, (k0 + 1) % 3);
        // std::cout << "    ip1" << std::endl;
        const int ip1 = _FUV(f1, (k1 + 1) % 3);
        // std::cout << "    iq1" << std::endl;
        const int iq1 = _FUV(f1, k1);
        // std::cout << "    ...." << std::endl;

        // uv coords of the edge as complex numbers
        const std::complex<double> p0(_UV(ip0, 0), _UV(ip0, 1));
        const std::complex<double> q0(_UV(iq0, 0), _UV(iq0, 1));
        const std::complex<double> p1(_UV(ip1, 0), _UV(ip1, 1));
        const std::complex<double> q1(_UV(iq1, 0), _UV(iq1, 1));

        // difference
        const std::complex<double> tp    = p1 - r01 * p0;
        const std::complex<double> tq    = q1 - r01 * q0;
        const std::complex<double> t_avg = 0.5 * tp + 0.5 * tq;

        double itu = std::round(t_avg.real());
        double itv = std::round(t_avg.imag());

        const double dtu = std::abs(itu - t_avg.real()); // in [0, 0.5]
        const double dtv = std::abs(itv - t_avg.imag()); // in [0, 0.5]

        if (itu == 1.0 && dtu > 0.49)
          itu = 0.;
        if (itv == 1.0 && dtv > 0.49)
          itv = 0.;

        // store
        _R_transfn(f0, k0)         = r01;
        _T_transfn(f0, k0)         = t_avg;
        _T_integer_transfn(f0, k0) = std::complex<double>(itu, itv);
        _T_decimal_transfn(f0, k0) = std::complex<double>(dtu, dtv);
        _flip_transfn(f0, k0)      = (pjump % 2) == 1;
      }
    }

    // debug : print the rounding error for constrained edges
    // for (int i = 0; i < ff__constrE_f0.rows(); i++) {
    //   const int f0     = ff__constrE_f0(i);
    //   const int k0     = ff__constrE_k0(i);
    //   const double dtu = _T_decimal_transfn(f0, k0).real();
    //   const double dtv = _T_decimal_transfn(f0, k0).imag();
    //   if (dtu > 1e-4 || dtv > 1e-4) { printf("f=%d, k=%d :  dtu= %0.4f, dtv= %0.4f\n", f0, k0, dtu, dtv); }
    // }

    return true;
  }

  // ===========================================================================

} // namespace IGSV
