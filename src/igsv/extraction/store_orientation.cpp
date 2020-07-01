#include <Eigen/Core>

#include <igsv/extraction/exact_predicates.h>
#include <igsv/extraction/store_orientation.h>

namespace IGSV {

  bool store_orientation(const Eigen::MatrixXd& _UV,                                //
                         const Eigen::MatrixXi& _FUV,                               //
                         const Eigen::Matrix<bool, Eigen::Dynamic, 1>& _F_isNarrow, //
                         Eigen::VectorXi& _F_uvOrientation) {

    _F_uvOrientation.resize(_FUV.rows(), 1);
    _F_uvOrientation.fill(1);

    for (int f = 0; f < _FUV.rows(); ++f)
      if (_F_isNarrow(f))
        _F_uvOrientation(f) =
            ORIENT2D(Point_2(_UV(_FUV(f, 0), 0), _UV(_FUV(f, 0), 1)), Point_2(_UV(_FUV(f, 1), 0), _UV(_FUV(f, 1), 1)),
                     Point_2(_UV(_FUV(f, 2), 0), _UV(_FUV(f, 2), 1)));

    return true;
  }

} // namespace IGSV
