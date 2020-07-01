//=============================================================================
//
//  STRUCT : Sample
//           ChainSample
//
//=============================================================================

#pragma once

//== NAMESPACES ===============================================================

namespace IGSV {

  //== STRUCT DEFINITION ====================================================

  struct Sample { // associated to a pixel

    Sample(const Sample& s)
        : pid(s.pid), qei(s.qei), t(s.t), d(s.d), px(s.px), py(s.py), tx(s.tx), ty(s.ty), ex(s.ex), ey(s.ey), w(s.w) {}

    Sample(int _pid, int _qei, double _t, double _d, double _px, double _py, double _tx, double _ty, double _ex,
           double _ey, double _w)
        : pid(_pid), qei(_qei), t(_t), d(_d), px(_px), py(_py), tx(_tx), ty(_ty), ex(_ex), ey(_ey), w(_w) {}

    // const properties
    const int pid;   // index of the pixel
    const int qei;   // index of a nearby q-edge (*not* necessarily the closest one)
    const double t;  // parameter of the closest point on the qedge
    const double d;  // distance to the closest point on the edge
    const double px; // pixel position x-coord
    const double py; // pixel position y-coord
    const double tx; // pixel tangent x-coord
    const double ty; // pixel tangent y-coord
    const double ex; // x-coord on the qedge: (1-t) * x0 + t * x1
    const double ey; // y-coord on the qedge: (1-t) * y0 + t * y1
    const double w;  // weight
  };

  //===========================================================================

  struct ChainSample : public Sample {

    ChainSample(const Sample& _sample, int _chain, double _t_chain)
        : Sample(_sample), chain(_chain), t_chain(_t_chain) {}

    const int chain;
    const double t_chain;

    // // fitting weight
    // double w() const { return pixel->w(); }
    //
    // // point on the curve for the parameter t_curve
    // void evaluate_bp();
    // std::complex<double> bp;
    //
    // // fitting cost, computed using w, mult, distance of pixel to bp
    // double cost() const;
  };

  //===========================================================================

  struct CurveSample : public Sample {

    CurveSample(const Sample& _sample, int _curve, double _t_curve)
        : Sample(_sample), curve(_curve), t_curve(_t_curve) {}

    const int curve;
    const double t_curve;
  };

  //===========================================================================

} // namespace IGSV
