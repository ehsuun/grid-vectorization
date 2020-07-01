#include <cmath>

#include <igsv/parametrization/pix_coords_from_xy.h>

namespace IGSV {
  void pix_coords_from_xy(double w, double h, double x, double y, int& i, int& j) {
    j = std::floor(x);
    i = h - 1 - std::floor(y);
    if (i < 0)
      i = 0;
    if (i > h - 1)
      i = h - 1;
    if (j < 0)
      j = 0;
    if (j > w - 1)
      j = w - 1;
  }
} // namespace IGSV
