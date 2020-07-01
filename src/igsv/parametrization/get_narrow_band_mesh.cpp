#include <Eigen/Core>

#include <igsv/parametrization/get_narrow_band_mesh.h>

#include <map>
#include <set>
#include <vector>

namespace IGSV {

  bool get_narrow_band_mesh(const Eigen::MatrixXd& _V,                                 //
                            const Eigen::MatrixXi& _F,                                 //
                            const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                            Eigen::MatrixXd& _V0,                                      //
                            Eigen::MatrixXi& _F0) {

    //// index maps between base mesh and narrow and sub-mesh
    std::map<int, int> V_base_to_sub, V_sub_to_base;
    std::map<int, int> F_base_to_sub, F_sub_to_base;
    std::map<int, int> U_base_to_sub, U_sub_to_base;
    std::set<int> v_used, u_used;

    // -- temp faces
    Eigen::MatrixXi F0_temp;
    // Eigen::MatrixXi FUV0_temp;
    {
      F0_temp.resize(_F.rows(), 3);
      // FUV0_temp.resize(_F.rows(), 3);

      int fs = 0;
      for (int fb = 0; fb < _F.rows(); ++fb)
        if (_F_isNarrow(fb)) {

          F0_temp.row(fs) << _F.row(fb);

          v_used.insert(_F(fb, 0));
          v_used.insert(_F(fb, 1));
          v_used.insert(_F(fb, 2));

          // u_used.insert(_FUV(fb, 0));
          // u_used.insert(_FUV(fb, 1));
          // u_used.insert(_FUV(fb, 2));

          F_base_to_sub.insert(std::make_pair(fb, fs));
          F_sub_to_base.insert(std::make_pair(fs, fb));

          fs++;
        }
      F0_temp.conservativeResize(fs, 3);
      // FUV0_temp.conservativeResize(fs, 3);
    }

    // -- vertices
    {
      _V0.resize(_V.rows(), 3);
      int vs = 0;
      for (auto vb : v_used) {
        V_base_to_sub.insert(std::make_pair(vb, vs));
        V_sub_to_base.insert(std::make_pair(vs, vb));
        _V0.row(vs++) = _V.row(vb);
      }
      _V0.conservativeResize(vs, 3);
    }

    // // -- uv vertices
    // {
    //   _UV0.resize(this->UV.rows(), 2);
    //   int us = 0;
    //   for (auto ub : u_used) {
    //     U_base_to_sub.insert(std::make_pair(ub, us));
    //     U_sub_to_base.insert(std::make_pair(us, ub));
    //     _UV0.row(us++) = this->UV.row(ub);
    //   }
    //   _UV0.conservativeResize(us, 2);
    // }

    // -- faces
    _F0.setZero(F0_temp.rows(), 3);
    // _FUV0.setZero(F0_temp.rows(), 3);
    for (int f = 0; f < F0_temp.rows(); ++f)
      for (int k = 0; k < 3; ++k) {
        _F0(f, k) = V_base_to_sub[F0_temp(f, k)];
        // _FUV0(f, k) = U_base_to_sub[FUV0_temp(f, k)];
      }

    return true;
  }

} // namespace IGSV
